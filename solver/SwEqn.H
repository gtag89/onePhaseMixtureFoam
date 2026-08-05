
Info << "SATURATION CORRECTION" << endl;

// Inverse of momentum equation diagonal coefficient 
volScalarField rAU(1.0/UEqn.A());

// Face-interpolated inverse diagonal coefficient
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

// Density-weighted inverse diagonal coefficient at faces
volScalarField rhoRAU("rhoRAU", rho*rAU);
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhoRAU));

// H/A term with capillary pressure constraint for momentum equation
volVectorField HbyA
(
    "HbyA", 
    constrainHbyA(rAU*UEqn.H(), U, pc)
);

// Gravitational flux contribution
surfaceScalarField phig
(
    "phig", 
    rhorAUf*fvc::interpolate(rho)*(g & mesh.Sf())
); 

// Calculate uncorrected mass flux using Rhie-Chow interpolation
// This prevents pressure-velocity decoupling on collocated grids
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        // Convective flux from uncorrected velocity
        fvc::flux(rho*HbyA)
        
        // Time derivative correction for momentum equation
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
      
        // Gravitational flux
      + phig
    )
);

// Capillary diffusion coefficient 
volScalarField diffGradSw
(
    "diffGradSw", 
    rho * rAU * dpcdSw * Sw
);

volScalarField Pe
(
    "Pe",
    (mag(U) * L/mesh.nCells()) / diffGradSw
); 

Info<< "Min. Peclet number = " << gMin(Pe) << " Max. Peclet number = " << gMax(Pe) << endl;
// Make flux relative to moving reference frame if present
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Face-interpolated diffusion coefficient for saturation equation
surfaceScalarField rhorAUfConstrain
(
    "rhorAUfConstrain", 
    rhorAUf*dpcdSwf
);

// Update saturation BCs to ensure flux consistency across domain
// This constrains the saturation field similar to pressure correction
constrainPressure(Sw, rho, U, phiHbyA, rhorAUfConstrain, MRF);

// Convert flux to relative form (remove mesh motion effects)
fvc::makeRelative(phiHbyA, rho, U);

// Construct saturation transport equation
fvScalarMatrix SwEqn
(
    // Temporal term: 
    fvm::ddt(rho, Sw) 
    
    // Convective term: 
  + fvc::div(phiHbyA, Sw)
  
    // Capillary diffusion term: 
  + fvm::laplacian(diffGradSw, Sw) 
  
    // Source terms (if any)
  ==
    fvOptions(rho, Sw)
);

while (pimple.correctNonOrthogonal())
{
    Info << "PERFORMING SATURATION CORRECTION" << endl;
    
    // Solve the saturation equation
    SwEqn.solve();
      
    // On final non-orthogonal iteration, update fluxes and velocity
    if (pimple.finalNonOrthogonalIter())
    {
        Info << "FLUX AND VELOCITY UPDATE" << endl;
        
        // Update conservative mass flux with saturation equation correction
        // This ensures mass conservation and pressure-velocity coupling
        phi = phiHbyA + SwEqn.flux()/fvc::interpolate(Sw);  
        
        // Explicitly relax saturation field for stability
        Sw.relax();
        
        // Correct velocity with pressure gradient from saturation equation
        // This couples capillary effects into the momentum balance
        U = HbyA + rAU*fvc::reconstruct((phig + SwEqn.flux()/fvc::interpolate(Sw))/rhorAUf);
        
        // Ensure velocity boundary conditions are satisfied
        U.correctBoundaryConditions();
        
        // Apply any finite volume source corrections
        fvOptions.correct(U);
    }
}

// Update porous media properties based on new saturation field
#include "updatePorousFields.H"

// Check saturation continuity and mass conservation
#include "saturationContinuityErrs.H"  

// Calculate inner loop convergence residual
scalar innerResidual = 
    (
        fvc::domainIntegrate(mag(Sw - SwPrev)) /
        fvc::domainIntegrate(SwPrev) 
    ).value();
    
// Store current saturation for next iteration
SwPrev = Sw;

Info<< "Inner loop saturation residual = " << innerResidual << endl;
if (innerResidual <= 1.0e-6)
{
    Info<< "****************   CONVERGED correction loop   ********************** " << endl;
}

// ************************************************************************* //
