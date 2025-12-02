
// Correct boundary velocities for moving reference frames
MRF.correctBoundaryVelocity(U);

Info << "ASSEMBLING MOMENTUM EQUATION" << endl;

fvVectorMatrix UEqn
( 
    fvm::ddt(rho, U) + fvm::div(phi, U) + MRF.DDt(rho, U)
    // Standard turbulence model (commented out for direct viscous modeling)
    // + turbulence->divDevRhoReff(U)
    // Laplacian viscous term: 
  - fvm::laplacian(mu, U)
    // Viscous stress divergence:  
    // This handles the non-isotropic part of the viscous stress tensor
  - fvc::div(mu * dev2(T(fvc::grad(U))))
    // Bulk viscosity gradient (explicit): 
  - fvc::grad(zeta*fvc::div(U))
    // Bulk viscosity Laplacian (implicit part):
  - fvm::laplacian(zeta, U)
    // Bulk viscosity Laplacian (explicit correction): 
  + fvc::laplacian(zeta, U)
    // Darcy resistance: 
  + fvm::Sp(dragTerm, U)
  ==
    // User-defined source terms (fvOptions)
    fvOptions(rho, U)
);

// Apply under-relaxation for stability in outer iterations
UEqn.relax();

// Apply any user-defined constraints to the equation
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    Info << "SOLVING MOMENTUM PREDICTOR WITH CAPILLARY PRESSURE" << endl;
    
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
              // Capillary pressure surface normal gradient: 
              - dpcdSwf*fvc::snGrad(Sw)
            ) * mesh.magSf()  // Convert to flux [kg/s]
            
            // Gravitational flux: 
          + fvc::interpolate(rho)*(g & mesh.Sf())
        )
    );
    
    // Apply any velocity corrections from fvOptions
    fvOptions.correct(U);
}

// ************************************************************************* //
