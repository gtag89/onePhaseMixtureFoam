/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

/*
Description
    Sets the reciprocal time-step (rDeltaT) for adaptive time stepping based 
    on local Courant number constraints for unsaturated porous media flow.
    
    This version is adapted from the standard compressible setRDeltaT.H with:
    - Compressibility effects (psi) removed for incompressible flow
    - Focus on convective time step limitations only
    
    The time step is controlled by:
    - Local Courant number (maxCo)
    - Maximum allowed time step (maxDeltaT)  
    - Smoothing to prevent rapid time step changes
    - Damping to limit rate of time step increase
*/

{
    volScalarField& rDeltaT = trDeltaT.ref();
    const dictionary& pimpleDict = pimple.dict();
    
    // Maximum Courant number for stability
    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.8)
    );
    
    // Smoothing coefficient for time step changes (0-1)
    // Lower values = more smoothing
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );
    
    // Damping coefficient for time step increase (0-1)  
    // Lower values = slower time step increases
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );
    
    // Maximum allowed time step [s]
    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );
    
    // Store current time step for damping calculations
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);
    
    // Calculate reciprocal time step based on local Courant number
    rDeltaT.ref() = max
    (
        // Reciprocal of maximum allowed time step
        1.0/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        
        // Convective time step limitation
        // Sum of face fluxes divided by (2*maxCo*cellVolume*density)
        fvc::surfaceSum(mag(phi))()()
        /((2.0*maxCo)*mesh.V()*rho())
    );
    
    /*
    // Transonic correction (disabled for incompressible flow)
    // This would add acoustic time step limitations for compressible flow
    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi)*fvc::flux(U)
        );
        rDeltaT.ref() = max
        (
            rDeltaT(),
            fvc::surfaceSum(mag(phid))()()
            /((2*maxCo)*mesh.V()*psi())
        );
    }
    */
    
    // Update boundary conditions for reciprocal time step field
    rDeltaT.correctBoundaryConditions();
    
    // Report initial time step range
    Info<< "Flow time scale min/max = "
        << gMin(1.0/rDeltaT.primitiveField()) << "/"
        << gMax(1.0/rDeltaT.primitiveField()) << " [s]" << endl;
    
    // Apply smoothing to reduce time step oscillations
    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
        
        Info<< "Smoothed flow time scale min/max = "
            << gMin(1.0/rDeltaT.primitiveField()) << "/"
            << gMax(1.0/rDeltaT.primitiveField()) << " [s]" << endl;
    }
    
    // Apply damping to limit rate of time step increase
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        // Allow unlimited decrease but limit increase rate
        rDeltaT = 
            rDeltaT0 * max
            (
                rDeltaT/rDeltaT0,                    // New/old ratio
                scalar(1.0) - rDeltaTDampingCoeff   // Maximum allowed decrease in rDeltaT
                                                     // (= maximum allowed increase in deltaT)
            );
            
        Info<< "Damped flow time scale min/max = "
            << gMin(1.0/rDeltaT.primitiveField()) << "/"
            << gMax(1.0/rDeltaT.primitiveField()) << " [s]" << endl;
    }
}

// ************************************************************************* //
