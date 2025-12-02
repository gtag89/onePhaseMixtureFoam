/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    onePhaseMixtureFoam

Group
    compressibleSolvers

Description
    Transient solver for unsaturated porous media flow using mixture theory.
    
    Solves coupled momentum and saturation transport equations for single-phase 
    water flow in variably saturated porous media with:
    
    * Van Genuchten capillary pressure relationships
    * Saturation-dependent relative permeability  
    * Darcy drag with porous media resistance
    * Bulk viscosity effects for numerical stability
    * Capillary pressure gradients as momentum source
    
    The solver uses a PIMPLE algorithm adapted for saturation-velocity coupling,
    where saturation gradients drive momentum transport through capillary forces.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "vanGenuchten.H"
#include "brooksCorey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for unsaturated porous media flow using mixture theory"
        " with van Genuchten capillary pressure relationships,"
        " adaptive time stepping and optional mesh motion."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createRhoUfIfPresent.H"
    
    // Note: Turbulence model validation disabled for direct viscous modeling
    // turbulence->validate();

    // Initialize time stepping for non-LTS cases
    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // =======================================================================
        // DYNAMIC MESH AND TIME STEP CONTROL
        // =======================================================================
        
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh for phi correction
        // Currently disabled - uncomment if mesh motion is needed
        /*
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU.reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }
        */

        // Set time step based on Courant number constraints
        if (LTS)
        {
            // Local time stepping - adaptive time step per cell
            #include "setRDeltaT.H"
        }
        else
        {
            // Global time stepping - single time step for entire domain
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Store velocity from previous time step for convergence monitoring
        volVectorField UPrev = U;
        
        // Outer PIMPLE loop for pressure-velocity coupling
        while (pimple.loop())
        {
            // Solve momentum equation with capillary pressure gradients
            #include "UEqn.H"
             
            // Store saturation from previous iteration for convergence monitoring
            volScalarField SwPrev = Sw;

            // Inner corrector loop for saturation-velocity coupling
            while (pimple.correct())
            {
                // Solve saturation transport with momentum coupling
                #include "SwEqn.H"
            }
        }
        
        // Update saturation history for next time step
        SwOld = Sw;

        // Write solution fields to disk
        runTime.write();

        // Report computational performance
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
