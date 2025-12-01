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
    {}
    
    // ==========================================================================
    // PUBLIC MEMBER FUNCTIONS
    // ==========================================================================
    
    //- Compute capillary pressure from effective saturation
    //  pc = pc0 * [(Se^(-1/m) - 1)^(1/n)]
    //  Units: [Pa]
    volScalarField computePc(const volScalarField& Se) const
    {
        return pc0_ * pow(pow(Se, -1.0/m_) - 1.0, 1.0/n_);
    }
    
    //- Compute derivative of capillary pressure with respect to water saturation
    //  dpc/dSw = -(pc0/(m*n)) * (1/(1-Sr)) * [(Se^(-1/m)-1)^(1/n-1)] * Se^(-1/m-1)
    //  Units: [Pa] (since dSw is dimensionless)
    volScalarField computedPcdSw(const volScalarField& Se, const dimensionedScalar Sr) const
    {
        return -(pc0_ / (m_ * n_)) * (1.0 / (1.0 - Sr)) * pow(pow(Se, -1.0/m_) 
               - 1.0, 1.0/n_ - 1.0) * pow(Se, -1.0/m_ - 1.0);
    }
    
    //- Compute relative permeability from effective saturation  
    //  Kr = Se^l * [1 - (1 - Se^(1/m))^m]^2
    //  This is the Mualem-van Genuchten model
    //  Units: [-] (dimensionless)
    volScalarField computeKr(const volScalarField& Se) const
    {
        return pow(Se, 0.5) * pow(1 - pow(1 - pow(Se, 1.0/m_), m_), 2);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
