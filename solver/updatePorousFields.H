/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

/*
Description
    Updates all porous media properties based on the current saturation field.
    
    This routine:
    1. Calculates effective saturation with numerical limits
    2. Updates hysteresis regime (currently disabled)
    3. Interpolates van Genuchten parameters for current regime
    4. Computes relative permeability and capillary pressure relationships
    5. Updates Darcy drag coefficient and bulk viscosity
    
    Called after each saturation equation solution to maintain consistency
    between saturation field and porous media properties.
*/

Info<< "------------- UPDATING POROUS MODEL VARIABLES -------------" << endl;

// =============================================================================
// EFFECTIVE SATURATION UPDATE
// =============================================================================

// Calculate effective saturation: Se = (Sw - Sr)/(Smax - Sr)
// For typical applications where Smax = 1: Se = (Sw - Sr)/(1 - Sr)
Se = (Sw - Sr)/(1.0 - Sr);

// Apply numerical bounds to prevent division by zero and ensure stability
// Limits Se to range [SeMin, 1-SeMin] for robust van Genuchten evaluation
Se = max(min(1.0 - SeMin, Se), SeMin);

// Update face-interpolated saturation for flux calculations
Swf = fvc::interpolate(Sw); 

// =============================================================================
// HYSTERESIS REGIME DETERMINATION (CURRENTLY DISABLED)
// =============================================================================

/*
// Hysteresis detection based on global saturation levels
// This would track whether the domain is in imbibition or drainage
if (gMax(Sw) < (Smax.value() - 0.01))
{
    // Domain saturation is increasing (imbibition phase)
    // Use imbibition parameters (mIm, nIm, pc0Im)
    Regime = 1.0;
}
else if (gMax(Sw) > (Smax.value() + 0.01))
{
    // Domain saturation is decreasing (drainage phase)  
    // Use drainage parameters (mDr, nDr, pc0Dr)
    Regime = 0.0;
}
else
{
    // Near equilibrium or mixed conditions
    // Default to imbibition parameters
    Regime = 1.0; 
}
*/

// =============================================================================
// VAN GENUCHTEN PARAMETER INTERPOLATION
// =============================================================================

// Interpolate van Genuchten parameters based on current hysteresis regime
// Regime = 1.0: pure imbibition, Regime = 0.0: pure drainage
m = Regime * mIm.value() + (1.0 - Regime) * mDr.value();
n = Regime * nIm.value() + (1.0 - Regime) * nDr.value();
l = Regime * lIm.value() + (1.0 - Regime) * lDr.value();
pc0 = Regime * pc0Im + (1.0 - Regime) * pc0Dr;

// =============================================================================
// VAN GENUCHTEN MODEL RECONSTRUCTION
// =============================================================================

// Create new van Genuchten model instance with updated parameters
vanGenuchten vg(m, n, l, pc0);

// =============================================================================
// RELATIVE PERMEABILITY UPDATE
// =============================================================================

// Compute relative permeability using van Genuchten model
// Kr = Se^l * [1 - (1 - Se^(1/m))^m]^2
Kr = vg.computeKr(Se) + KrSmall;

// KrSmall provides numerical stability for very low permeabilities
// Prevents division by zero in Darcy drag calculations

// =============================================================================
// CAPILLARY PRESSURE RELATIONSHIPS
// =============================================================================

// Compute capillary pressure: pc = pc0 * [(Se^(-1/m) - 1)^(1/n)]
pc = vg.computePc(Se); 

// Compute derivative dpc/dSw for diffusion coefficient
// This is needed for the capillary diffusion term in saturation equation
dpcdSw = vg.computedPcdSw(Se, Sr);

// Apply numerical limiters to prevent excessive gradients
dpcdSw = max(min(dpcdSw, dpcdSwMax), dpcdSwMin);

// Interpolate derivative to faces for flux calculations
dpcdSwf = fvc::interpolate(dpcdSw); 

// =============================================================================
// FLOW RESISTANCE UPDATES
// =============================================================================

// Update Darcy drag coefficient [kg/(m³·s)]
// dragTerm = (ε·Sw·μ)/(Kr·Kabs)
// This represents the momentum loss due to viscous flow through porous media
dragTerm = (eps * Sw * mu)/(Kr * Kabs);

// =============================================================================
// BULK VISCOSITY UPDATE  
// =============================================================================

// Reset bulk viscosity to base coefficient
// This could be modified to depend on saturation or other properties
zeta = zetaCoeff;

// Alternative saturation-dependent bulk viscosity models (examples):
// zeta = zetaCoeff * pow(1.0 - Sw, alpha);                    // Power law
// zeta = zetaCoeff * (1.0 + alpha * mag(fvc::grad(Sw)));      // Gradient-dependent  
// zeta = zetaCoeff * pow(1.0 - Sw, alpha) * pow(Sw, beta);    // Double power law

Info<< "Porous media properties updated successfully" << endl;

// ************************************************************************* //

Info << "Saturation: " << " Min(S) = " << gMin(Sw) << " Max(S) = " << gMax(Sw) << endl;
//Info << "Saturation upstream: " << Sw[3197] << endl;
Info << "Capillary pressure : " << " Min(pc) = " << gMin(pc) << " Max(pc) = " << gMax(pc) << endl;
Info << "dpcdSw: " << " Min(dpcdSw) = " << gMin(dpcdSw) << " Max(dpcdSw) = " << gMax(dpcdSw) << endl;
Info << "Relative permeability: " << " Min(Kr) = " << gMin(Kr) << " Max(Kr) = " << gMax(Kr) << endl;
Info << "Sec Viscosity: " << " Min(zeta) = " << gMin(zeta) << " Max(zeta) = " << gMax(zeta) << endl; 
Info << "Velocity from UEqn: " << " Min(U) = " << gMin(U) << " Max(U) = " << gMax(U) << endl;

