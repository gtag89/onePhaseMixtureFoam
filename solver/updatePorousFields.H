
Info<< "------------- UPDATING POROUS MODEL VARIABLES -------------" << endl;

// Calculate effective saturation: Se = (Sw - Sr)/(Smax - Sr)
// For typical applications where Smax = 1: Se = (Sw - Sr)/(1 - Sr)
Se = (Sw - Sr)/(1.0 - Sr);

// Apply numerical bounds to prevent division by zero and ensure stability
// Limits Se to range [SeMin, 1-SeMin] for robust van Genuchten evaluation
Se = max(min(1.0 - SeMin, Se), SeMin);

// Update face-interpolated saturation for flux calculations
Swf = fvc::interpolate(Sw); 

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

// Interpolated parameters based on current regime
dimensionedScalar n = Regime * nIm.value() + (1.0 - Regime) * nDr.value();
dimensionedScalar l = Regime * lIm.value() + (1.0 - Regime) * lDr.value();
volScalarField pc0 = Regime * pc0Im + (1.0 - Regime) * pc0Dr;
dimensionedScalar m = Regime * mIm.value() + (1.0 - Regime) * mDr.value();
dimensionedScalar lambda = Regime * lambdaIm.value() + (1.0 - Regime) * lambdaDr.value();
volScalarField pcd = Regime * pcdIm + (1.0 - Regime) * pcdDr;
dimensionedScalar eta = Regime * etaIm.value() + (1.0 - Regime) * etaDr.value();

// Porous quantities update
if (porousModel == "vanGenuchten")
{
    
    vanGenuchten vg(m, n, l, pc0);
    
    Kr = vg.computeKr(Se) + KrSmall;
    pc = vg.computePc(Se);
    dpcdSw = vg.computedPcdSw(Se, Sr);
}
else if (porousModel == "brooksCorey")
{

    brooksCorey bc(lambda, eta, pcd);
    
    Kr = bc.computeKr(Se) + KrSmall;
    pc = bc.computePc(Se);
    dpcdSw = bc.computedPcdSw(Se, Sr);
}
else // porousModel == "none"
{
    Info<< "No porous model - using default values" << nl << endl;
}

// Apply numerical limiters to prevent excessive gradients
dpcdSw = max(min(dpcdSw, dpcdSwMax), dpcdSwMin);

// Interpolate derivative to faces for flux calculations
dpcdSwf = fvc::interpolate(dpcdSw); 

// Update Darcy drag coefficient [kg/(m³·s)]
dragTerm = (eps * Sw * mu)/(Kr * Kabs);

// Reset bulk viscosity to base coefficient
// This could be modified to depend on saturation or other properties
//zeta = zetaCoeff;

// Alternative saturation-dependent bulk viscosity models (examples):
zeta = zetaCoeff * (1.0 + alpha * exp(Sw/(beta * mag(1.0 - Sw))));

Info<< "Porous media properties updated successfully" << endl;

Info << "Saturation: " << " Min(S) = " << gMin(Sw) << " Max(S) = " << gMax(Sw) << endl;

Info << "Capillary pressure : " << " Min(pc) = " << gMin(pc) << " Max(pc) = " << gMax(pc) << endl;
Info << "dpcdSw: " << " Min(dpcdSw) = " << gMin(dpcdSw) << " Max(dpcdSw) = " << gMax(dpcdSw) << endl;
Info << "Relative permeability: " << " Min(Kr) = " << gMin(Kr) << " Max(Kr) = " << gMax(Kr) << endl;
Info << "Sec Viscosity: " << " Min(zeta) = " << gMin(zeta) << " Max(zeta) = " << gMax(zeta) << endl; 
Info << "Velocity from UEqn: " << " Min(U) = " << gMin(U) << " Max(U) = " << gMax(U) << endl;

