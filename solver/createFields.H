/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

// Create time-step control for adjustable time stepping
#include "createRDeltaT.H"

// =============================================================================
// VELOCITY FIELD
// =============================================================================
Info<< "Reading velocity field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// =============================================================================
// THERMOPHYSICAL AND TRANSPORT PROPERTIES
// =============================================================================
Info << nl << "Reading thermophysical properties" << endl;
IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Extract mixture and transport properties
dictionary mixture(thermophysicalProperties.subDict("mixture"));
dictionary transport(mixture.subDict("transport"));

// Water viscosity [kg/(m·s)]
dimensionedScalar muw
(    
    "muw",
    dimensionSet(1, -1, -1, 0, 0, 0, 0), 
    transport.lookupOrDefault("muw", 0.001)  // Default to water viscosity at 20°C
);

// Water density [kg/m³]
dimensionedScalar rhow
(    
    "rhow",
    dimensionSet(1, -3, 0, 0, 0, 0, 0), 
    transport.lookupOrDefault("rhow", 1000.0)  // Default to water density
);

// Characteristic length scale [m]
dimensionedScalar L
(    
    "L",
    dimensionSet(0, 1, 0, 0, 0, 0, 0), 
    transport.lookupOrDefault("L", 1.0)
);

// Zeta coefficient for non-Darcian flow [kg/(m·s)]
dimensionedScalar zetaCoeff
(    
    "zetaCoeff",
    dimensionSet(1, -1, -1, 0, 0, 0, 0),
    transport.lookupOrDefault("zetaCoeff", 0.0)
);

// Alpha parameter (dimensionless) - model dependent
dimensionedScalar alpha
(    
    "alpha",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),  
    transport.lookupOrDefault("alpha", 1.0)
);

// Beta parameter (dimensionless) - model dependent
dimensionedScalar beta
(    
    "beta",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),  
    transport.lookupOrDefault("beta", 1.0)
);

// =============================================================================
// GRAVITY AND REFERENCE FIELDS
// =============================================================================
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// =============================================================================
// SATURATION FIELDS
// =============================================================================
Info<< "Reading water saturation field Sw\n" << endl;
volScalarField Sw
(
    IOobject
    (
        "Sw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Store previous time step saturation for time derivatives
volScalarField SwOld = Sw;

// =============================================================================
// COMPUTED FIELDS
// =============================================================================

// Dynamic viscosity field [kg/(m·s)]
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    muw
);

// Include porous media specific fields
#include "createPorousFields.H"

// Density field [kg/m³]
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rhow
);

// Create mass flux field
#include "compressibleCreatePhi.H"

// =============================================================================
// NON-DARCIAN FLOW COEFFICIENT
// =============================================================================
const label& nCells = mesh.nCells();

volScalarField zeta
(
    IOobject
    (
        "zeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeta", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
);

// Initialize zeta coefficient (various models available - see comments below)
zeta = zetaCoeff;

// Alternative zeta models (commented out - uncomment as needed):
// Forchheimer-type models:
// zeta = zetaCoeff * pow(1.0 - Sw, alpha);                                    // Power law with (1-Sw)
// zeta = zetaCoeff * pow(1.0 - Sw, alpha) * pow(Sw, beta);                   // Double power law
// 
// Capillary pressure dependent models:
// zeta = zetaCoeff * (1.0 + alpha * mag(dpcdSw/pc0));                        // Linear pc dependence
// zeta = zetaCoeff * (1.0 + alpha * pow(1.0/mag(dpcdSw/pc0), beta));         // Power law pc dependence
//
// Exponential models:
// zeta = zetaCoeff * (1.0 + alpha * exp(-pow((Sw - Smax)/Sw, beta)));        // Exponential transition
//
// Gradient dependent models:
// zeta = zetaCoeff * (1.0 + alpha * mag(fvc::grad(Sw)));                     // Saturation gradient dependent

/*
// Alternative: Cell-by-cell zeta assignment based on saturation thresholds
forAll(Sw, celli)
{
    if (Sw[celli] < alpha.value() || Sw[celli] > beta.value())
    {
        zeta[celli] = zetaCoeff.value();
    }
    else
    {
        zeta[celli] = zetaCoeff.value() * (1.0 + gamma.value() * Sw[celli]);
    }
}
*/

// =============================================================================
// SOLVER SETUP
// =============================================================================
// Set flux requirement for saturation equation
mesh.setFluxRequired(Sw.name());

// Include MRF (Multiple Reference Frame) zones if present
#include "createMRF.H"

// Include finite volume options (fvOptions) for source terms
#include "createFvOptions.H"

