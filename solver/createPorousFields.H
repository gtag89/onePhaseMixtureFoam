
Info<< "------------- CREATING POROUS MODEL VARIABLES -------------" << endl;

Info<< "Reading porous models dictionary\n" << endl;
IOdictionary porousProperties
(
    IOobject
    (
        "porousProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading absolute permeability field Kabs [m²]\n" << endl;
volScalarField Kabs
(
    IOobject
    (
        "Kabs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading porosity field eps [-]\n" << endl;
volScalarField eps
(   
    IOobject
    (   
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

word porousModel(porousProperties.lookupOrDefault<word>("porousModel", "none"));

if (porousModel != "vanGenuchten" && porousModel != "brooksCorey" 
    && porousModel != "none")
{
    FatalErrorInFunction
        << "Unknown porous model: " << porousModel << nl
        << "Valid options are: vanGenuchten, brooksCorey, none" << nl
        << exit(FatalError);
}

Info<< "Selected porous model: " << porousModel << nl << endl;

// Handle case when no porous model is specified
if (porousModel == "none")
{
    Info<< "No porous model specified - setting porosity to unity\n" << endl;
    forAll(eps, celli)
    {
        eps[celli] = 1.0;
    }
}

// VAN GENUCHTEN MODEL COEFFICIENTS
dictionary coefficients(porousProperties.subDict("coefficients"));

// Residual saturation [-]
dimensionedScalar Sr
(
    "Sr", 
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("Sr", 0.0)
);

// Maximum saturation [-]
dimensionedScalar Smax
(
    "Smax", 
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("Smax", 1.0)
);

// Minimum effective saturation for numerical stability [-]
dimensionedScalar SeMin
(
    "SeMin", 
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("SeMin", 1e-6)
);

// Van Genuchten parameter n for imbibition [-]
dimensionedScalar nIm
(    
    "nIm",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("nIm", 2.0)
);

// Pore connectivity parameter for imbibition [-]
dimensionedScalar lIm
(    
    "lIm",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("lIm", 0.5)
);

// Air-entry pressure for imbibition [Pa]
volScalarField pc0Im
(   
    IOobject
    (   
        "pc0Im",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "pc0Im",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        coefficients.lookupOrDefault("pc0Im", 1000.0)
    )
);

// Van Genuchten parameter n for drainage [-]
dimensionedScalar nDr
(    
    "nDr",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("nDr", 2.0)
);

// Pore connectivity parameter for drainage [-]
dimensionedScalar lDr
(    
    "lDr",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("lDr", 0.5)
);

// Air-entry pressure for drainage [Pa]
volScalarField pc0Dr
(   
    IOobject
    (   
        "pc0Dr",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "pc0Dr",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        coefficients.lookupOrDefault("pc0Dr", 1000.0)
    )
);

dimensionedScalar mIm = 1.0 - 1.0/nIm; 
dimensionedScalar mDr = 1.0 - 1.0/nDr; 

// BROOKS-COREY MODEL COEFFICIENTS
// Brooks-Corey parameter lambda for imbibition [-]
dimensionedScalar lambdaIm
(    
    "lambdaIm",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("lambdaIm", 2.0)
);
// Air-entry pressure for imbibition [Pa]
volScalarField pcdIm
(   
    IOobject
    (   
        "pcdIm",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "pcdIm",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        coefficients.lookupOrDefault("pcdIm", 1000.0)
    )
);

// Brooks-Corey parameter lambda for drainage [-]
dimensionedScalar lambdaDr
(    
    "lambdaDr",
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    coefficients.lookupOrDefault("lambdaDr", 2.0)
);
// Air-entry pressure for drainage [Pa]
volScalarField pcdDr
(   
    IOobject
    (   
        "pcdDr",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "pcdDr",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        coefficients.lookupOrDefault("pcdDr", 1000.0)
    )
);

dimensionedScalar etaIm = (2.0 + 3.0*lambdaIm)/lambdaIm; 
dimensionedScalar etaDr = (2.0 + 3.0*lambdaDr)/lambdaDr; 

// NUMERICAL LIMITERS
dictionary limiters(porousProperties.subDict("limiters"));

// Maximum allowed dpcdSw [Pa] (most negative value)
dimensionedScalar dpcdSwMax
(
    "dpcdSwMax", 
    dimensionSet(1, -1, -2, 0, 0, 0, 0), 
    limiters.lookupOrDefault("dpcdSwMax", -1e-6)
);

// Minimum allowed dpcdSw [Pa] (most positive value) 
dimensionedScalar dpcdSwMin
(
    "dpcdSwMin", 
    dimensionSet(1, -1, -2, 0, 0, 0, 0), 
    limiters.lookupOrDefault("dpcdSwMin", -1e+6)
);

// Minimum relative permeability for numerical stability [-]
dimensionedScalar KrSmall
(
    "KrSmall", 
    dimensionSet(0, 0, 0, 0, 0, 0, 0), 
    limiters.lookupOrDefault("KrSmall", 1e-6)
);

// HYSTERESIS MODEL
// Regime scalar: 1.0 = imbibition, 0.0 = drainage
// This should be computed based on saturation history
scalar Regime = 1.0;

// Interpolated parameters based on current regime
dimensionedScalar n = Regime * nIm.value() + (1.0 - Regime) * nDr.value();
dimensionedScalar l = Regime * lIm.value() + (1.0 - Regime) * lDr.value();
volScalarField pc0 = Regime * pc0Im + (1.0 - Regime) * pc0Dr;
dimensionedScalar m = Regime * mIm.value() + (1.0 - Regime) * mDr.value();
dimensionedScalar lambda = Regime * lambdaIm.value() + (1.0 - Regime) * lambdaDr.value();
volScalarField pcd = Regime * pcdIm + (1.0 - Regime) * pcdDr;
dimensionedScalar eta = Regime * etaIm.value() + (1.0 - Regime) * etaDr.value();

// Effective saturation [-]
volScalarField Se
(
    IOobject
    (
        "Se",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Sw - Sr)/(1.0 - Sr)
);

// Apply numerical limits to effective saturation
Se = max(min(1.0 - SeMin, Se), SeMin);

// Interpolated saturation at faces
surfaceScalarField Swf("Swf", fvc::interpolate(Sw));

// Relative permeability [-]
Info<< "Creating relative permeability field Kr [-]\n" << endl;
volScalarField Kr
(
    IOobject
    (
        "Kr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Kr", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
);

// Capillary pressure [Pa]
Info<< "Creating capillary pressure field pc [Pa]\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("pc", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0)
);

// Capillary pressure derivative [Pa]
Info<< "Creating dpcdSw field [Pa]\n" << endl;
volScalarField dpcdSw
(
    IOobject
    (
        "dpcdSw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dpcdSw", dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0)
);

if (porousModel == "vanGenuchten")
{
    Info<< "Initializing van Genuchten model" << nl << endl;

    vanGenuchten vg(m, n, l, pc0);
    
    Kr = vg.computeKr(Se) + KrSmall;
    pc = vg.computePc(Se);
    dpcdSw = vg.computedPcdSw(Se, Sr);
}
else if (porousModel == "brooksCorey")
{
    Info<< "Initializing Brooks-Corey model" << nl << endl;

    brooksCorey bc(lambda, eta, pcd);
    
    Kr = bc.computeKr(Se) + KrSmall;
    pc = bc.computePc(Se);
    dpcdSw = bc.computedPcdSw(Se, Sr);
}
else // porousModel == "none"
{
    Info<< "No porous model - using default values" << nl << endl;
}

// Apply numerical limiters
dpcdSw = max(min(dpcdSw, dpcdSwMax), dpcdSwMin);

// Interpolated derivative at faces
surfaceScalarField dpcdSwf("dpcdSwf", fvc::interpolate(dpcdSw)); 

// Darcy drag coefficient [kg/(m³·s)]
Info<< "Creating Darcy drag term [kg/(m³·s)]\n" << endl;
volScalarField dragTerm
(
    IOobject
    (
        "dragTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (eps * Sw * mu) / (Kr * Kabs)
);

