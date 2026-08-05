// Create time-step control for adjustable time stepping
#include "createRDeltaT.H"

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
const label inletPatchID = mesh.boundaryMesh().findPatchID("top");
if (inletPatchID < 0)
{
    FatalError << "Patch 'top' not found" << exit(FatalError);
}
// Face-area-weighted average of pore velocity magnitude at inlet
scalar UinVal = gAverage(mag(U.boundaryField()[inletPatchID]));
dimensionedScalar UinDim("UinDim", dimVelocity, UinVal);

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
    transport.lookupOrDefault("alpha", 0.0)
);

// Beta parameter (dimensionless) - model dependent
dimensionedScalar beta
(    
    "beta",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),  
    transport.lookupOrDefault("beta", 0.0)
);

// Gamma parameter (dimensionless) - model dependent
dimensionedScalar gamma
(    
    "gamma",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),  
    transport.lookupOrDefault("gamma", 0.0)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

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

// Include porous media specific fields
#include "createPorousFields.H"

// Bulk (second) viscosity field [kg/(m·s)]
volScalarField zeta
(
    IOobject
    (
        "zeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeta", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
);

// Initialize the bulk viscosity model from the inlet saturation and velocity
dimensionedScalar Sin("Sin", dimless, gAverage(Sw.boundaryField()[inletPatchID]));
dimensionedScalar avgU("avgU", dimless, UinVal);
dimensionedScalar zetaFunc = zetaCoeff * (1.0 + 70.3035 * pow(1.0 - Sin, 4.9699));
dimensionedScalar alphaFunc = alpha / (1.0 + exp(3.1143 * (log10(avgU) + 3.1213)));
volScalarField customFunc = exp(alphaFunc * 0.5 * Se / max(1 - Se, Sr));
zeta = zetaFunc * (1 + alphaFunc * customFunc);

// Set flux requirement for saturation equation
mesh.setFluxRequired(Sw.name());

// Include MRF (Multiple Reference Frame) zones if present
#include "createMRF.H"

// Include finite volume options (fvOptions) for source terms
#include "createFvOptions.H"

