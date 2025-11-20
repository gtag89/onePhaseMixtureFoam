/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef vanGenuchten_H
#define vanGenuchten_H

#include "volFields.H"
#include "dimensionedScalar.H"

/*---------------------------------------------------------------------------*\
                       Class vanGenuchten Declaration
\*---------------------------------------------------------------------------*/

/*
Description
    Van Genuchten model for unsaturated porous media relationships.
    
    Implements the classic van Genuchten (1980) model for:
    - Soil water retention curve (capillary pressure vs saturation)
    - Relative permeability as function of effective saturation
    - Derivative relationships for numerical solvers
    
    Model equations:
    - Capillary pressure: pc = pc0 * [(Se^(-1/m) - 1)^(1/n)]
    - Relative permeability: Kr = Se^l * [1 - (1 - Se^(1/m))^m]^2
    - Derivative: dpc/dSw = -(pc0/(m*n)) * (1/(1-Sr)) * [(Se^(-1/m)-1)^(1/n-1)] * Se^(-1/m-1)
    
    where:
    - Se = effective saturation = (Sw - Sr)/(Smax - Sr)
    - m, n = van Genuchten shape parameters (m = 1 - 1/n typically)
    - l = pore connectivity parameter (typically 0.5)
    - pc0 = air-entry pressure
    - Sr = residual saturation

References:
    van Genuchten, M.Th. (1980). A closed-form equation for predicting the 
    hydraulic conductivity of unsaturated soils. Soil Science Society of 
    America Journal, 44(5), 892-898.
*/

class vanGenuchten
{
private:
    
    // ==========================================================================
    // PRIVATE DATA MEMBERS
    // ==========================================================================
    
    //- Van Genuchten shape parameter m [-]
    //  Controls the steepness of the retention curve
    //  Typical range: 0.1 - 0.9
    dimensionedScalar m_;
    
    //- Van Genuchten shape parameter n [-]
    //  Related to pore size distribution
    //  Constraint: n > 1, typically n = 1/(1-m)
    dimensionedScalar n_;
    
    //- Pore connectivity parameter l [-]
    //  Controls relative permeability curve shape
    //  Typical value: 0.5 (Mualem assumption)
    dimensionedScalar l_;
    
    //- Air-entry pressure field pc0 [Pa]
    //  Pressure required to desaturate largest pores
    //  Can be spatially variable for heterogeneous media
    volScalarField pc0_;

public:
    
    // ==========================================================================
    // CONSTRUCTORS
    // ==========================================================================
    
    //- Constructor from components
    vanGenuchten
    (
        const dimensionedScalar m,
        const dimensionedScalar n, 
        const dimensionedScalar l,
        const volScalarField& pc0
    )
    :
        m_(m), 
        n_(n), 
        l_(l), 
        pc0_(pc0)
    {
        // Basic parameter validation
        if (m_.value() <= 0.0 || m_.value() >= 1.0)
        {
            FatalErrorInFunction
                << "van Genuchten parameter m must be in range (0,1)" << nl
                << "Provided value: " << m_.value() << nl
                << exit(FatalError);
        }
        
        if (n_.value() <= 1.0)
        {
            FatalErrorInFunction
                << "van Genuchten parameter n must be > 1" << nl
                << "Provided value: " << n_.value() << nl
                << exit(FatalError);
        }
        
        if (l_.value() < 0.0)
        {
            FatalErrorInFunction
                << "Pore connectivity parameter l must be >= 0" << nl
                << "Provided value: " << l_.value() << nl
                << exit(FatalError);
        }
    }
    
    // ==========================================================================
    // PUBLIC MEMBER FUNCTIONS
    // ==========================================================================
    
    //- Compute capillary pressure from effective saturation
    //  pc = pc0 * [(Se^(-1/m) - 1)^(1/n)]
    //  Units: [Pa]
    volScalarField computePc(const volScalarField& Se) const
    {
        // Ensure Se is within valid bounds to prevent numerical issues
        volScalarField SeBounded = max(min(Se, 0.9999), 1e-6);
        
        return pc0_ * pow
        (
            pow(SeBounded, -1.0/m_) - 1.0, 
            1.0/n_
        );
    }
    
    //- Compute derivative of capillary pressure with respect to water saturation
    //  dpc/dSw = -(pc0/(m*n)) * (1/(1-Sr)) * [(Se^(-1/m)-1)^(1/n-1)] * Se^(-1/m-1)
    //  Units: [Pa] (since dSw is dimensionless)
    volScalarField computedPcdSw
    (
        const volScalarField& Se, 
        const dimensionedScalar Sr
    ) const
    {
        // Ensure Se is within valid bounds
        volScalarField SeBounded = max(min(Se, 0.9999), 1e-6);
        
        // Pre-calculate common terms for efficiency
        volScalarField SeToMinusInvM = pow(SeBounded, -1.0/m_);
        volScalarField bracket = SeToMinusInvM - 1.0;
        
        return -(pc0_ / (m_ * n_)) 
             * (1.0 / (1.0 - Sr)) 
             * pow(bracket, 1.0/n_ - 1.0) 
             * pow(SeBounded, -1.0/m_ - 1.0);
    }
    
    //- Compute relative permeability from effective saturation  
    //  Kr = Se^l * [1 - (1 - Se^(1/m))^m]^2
    //  This is the Mualem-van Genuchten model
    //  Units: [-] (dimensionless)
    volScalarField computeKr(const volScalarField& Se) const
    {
        // Ensure Se is within valid bounds
        volScalarField SeBounded = max(min(Se, 0.9999), 1e-6);
        
        // Pre-calculate terms for numerical stability
        volScalarField SeToInvM = pow(SeBounded, 1.0/m_);
        volScalarField innerTerm = 1.0 - SeToInvM;
        volScalarField bracketed = 1.0 - pow(innerTerm, m_);
        
        return pow(SeBounded, l_) * pow(bracketed, 2.0);
    }
    
    // ==========================================================================
    // ACCESSOR FUNCTIONS
    // ==========================================================================
    
    //- Return van Genuchten parameter m
    const dimensionedScalar& m() const { return m_; }
    
    //- Return van Genuchten parameter n  
    const dimensionedScalar& n() const { return n_; }
    
    //- Return pore connectivity parameter l
    const dimensionedScalar& l() const { return l_; }
    
    //- Return air-entry pressure field
    const volScalarField& pc0() const { return pc0_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
