    
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

/*
Description
    Momentum equation for unsaturated porous media flow with:
    - Dynamic viscosity (first viscosity μ) 
    - Bulk viscosity (second viscosity ζ) with implicit-explicit corrections
    - Darcy drag in porous media
    - Capillary pressure gradients as driving force
    - Gravity effects
    
    The momentum equation is:
    ∂(ρU)/∂t + ∇·(ρUU) = -∇pc + ∇·τ + ρg - dragTerm·U + Sources
    
    where:
    - τ = μ(∇U + ∇U^T) + ζ∇·U·I (stress tensor with bulk viscosity)
    - ζ is the second viscosity coefficient (bulk viscosity)
    - dragTerm = εSμ/(KrKabs) (Darcy resistance)
    - pc = capillary pressure
*/

// =============================================================================
// BOUNDARY CONDITION PREPARATION
// =============================================================================

// Correct boundary velocities for moving reference frames
MRF.correctBoundaryVelocity(U);

Info << "ASSEMBLING MOMENTUM EQUATION" << endl;

// =============================================================================
// MOMENTUM EQUATION ASSEMBLY
// =============================================================================

fvVectorMatrix UEqn
(
    // ==========================================================================
    // TEMPORAL AND CONVECTIVE TERMS
    // ==========================================================================
    
    // Temporal term: ∂(ρU)/∂t
    fvm::ddt(rho, U) 
    
    // Convective term: ∇·(ρUU)  
  + fvm::div(phi, U) 
  
    // Moving reference frame correction
  + MRF.DDt(rho, U)
  
    // ==========================================================================
    // TURBULENCE MODEL (DISABLED)
    // ==========================================================================
    
    // Standard turbulence model (commented out for direct viscous modeling)
    // + turbulence->divDevRhoReff(U)
    
    // ==========================================================================
    // DYNAMIC VISCOSITY (FIRST VISCOSITY) EFFECTS
    // ==========================================================================
    
    // Laplacian viscous term: -∇·(μ∇U)
  - fvm::laplacian(mu, U) 
  
    // Viscous stress divergence: -∇·(μ∇U^T) 
    // This handles the non-isotropic part of the viscous stress tensor
  - fvc::div(mu * dev2(T(fvc::grad(U))))
  
    // ==========================================================================
    // BULK VISCOSITY (SECOND VISCOSITY) WITH IMPLICIT-EXPLICIT CORRECTIONS
    // ==========================================================================
    
    // Bulk viscosity gradient (explicit): -∇(ζ∇·U)
    // Represents resistance to compression/dilation
  - fvc::grad(zeta*fvc::div(U))
  
    // Bulk viscosity Laplacian (implicit part): -∇·(ζ∇U)
    // Treated implicitly for diagonal dominance and stability
  - fvm::laplacian(zeta, U)
  
    // Bulk viscosity Laplacian (explicit correction): +∇·(ζ∇U) 
    // Explicit correction to maintain physical accuracy while preserving stability
    // This implements: -∇(ζ∇·U) - ∇·(ζ∇U) + ∇·(ζ∇U) = -∇(ζ∇·U)
  + fvc::laplacian(zeta, U)
  
    // ==========================================================================
    // POROUS MEDIA RESISTANCE (DARCY DRAG)
    // ==========================================================================
    
    // Darcy resistance: dragTerm·U where dragTerm = εSμ/(KrKabs)
    // This represents viscous resistance in porous media
  + fvm::Sp(dragTerm, U)
  
    // ==========================================================================
    // SOURCE TERMS
    // ==========================================================================
    
  ==
    // User-defined source terms (fvOptions)
    fvOptions(rho, U)
);

// =============================================================================
// EQUATION CONDITIONING
// =============================================================================

// Apply under-relaxation for stability in outer iterations
UEqn.relax();

// Apply any user-defined constraints to the equation
fvOptions.constrain(UEqn);

// =============================================================================
// MOMENTUM PREDICTOR SOLUTION
// =============================================================================

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
                // =============================================================
                // CAPILLARY PRESSURE GRADIENT FORCE
                // =============================================================
                
                // Capillary pressure surface normal gradient: -∇pc·n̂
                // This is the primary driving force in unsaturated flow
                // dpcdSwf = ∂pc/∂Sw at faces, fvc::snGrad(Sw) = ∇Sw·n̂ at faces
                - dpcdSwf*fvc::snGrad(Sw)
                
                // =============================================================
                // GRAVITATIONAL FORCE  
                // =============================================================
                
            ) * mesh.magSf()  // Convert to flux [kg/s]
            
            // Gravitational flux: ρg·n̂·A
          + fvc::interpolate(rho)*(g & mesh.Sf())
        )
    );
    
    // Apply any velocity corrections from fvOptions
    fvOptions.correct(U);
}

// ************************************************************************* //