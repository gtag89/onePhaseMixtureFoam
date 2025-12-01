/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef brooksCorey_H
#define brooksCorey_H

#include "volFields.H"
#include "dimensionedScalar.H"

/*---------------------------------------------------------------------------*\
                       Class brooksCorey Declaration
\*---------------------------------------------------------------------------*/

/*
Description
    Brooks-Corey model for unsaturated porous media relationships.
    
    Model equations:
    - Capillary pressure: pc = pc0 * Se^(-1/n)
    - Derivative: dpc/dSw = -pc0 / (n * (1 - Sr)) * Se^(-1/n - 1)
    - Relative permeability: Kr = Se^m
      where m = (2 + 3*n)/n for Burdine's theory
    
    where:
    - Se = effective saturation = (Sw - Sr)/(1 - Sr)
    - n = pore size distribution index (typically 0.2 - 5.0)
    - pc0 = displacement (air-entry) pressure [Pa]
    - Sr = residual saturation
    - Smax = maximum saturation

References:
    Brooks, R.H. and Corey, A.T. (1964). Hydraulic properties of porous media.
    Hydrology Papers, Colorado State University, Fort Collins, Colorado.
*/

class brooksCorey
{
private:
    
    // ==========================================================================
    // PRIVATE DATA MEMBERS
    // ==========================================================================
    
    
    //- Pore size distribution index lambda [-]
    dimensionedScalar n_;
    
    //- Relative permeability exponent eta [-]
    dimensionedScalar m_;
    
    //- Displacement (air-entry) pressure field pd [Pa]
    volScalarField pc0_;

public:
    
    // ==========================================================================
    // CONSTRUCTORS
    // ==========================================================================

    brooksCorey
    (
        const dimensionedScalar n,
        const dimensionedScalar m, 
        const volScalarField& pc0
    )
    :
        n_(n), 
        m_(m), 
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
        return pc0_ * pow(Se, -1.0/n_);
    }
    
    //- Compute derivative of capillary pressure with respect to water saturation
    //  dpc/dSw = -(pc0/(m*n)) * (1/(1-Sr)) * [(Se^(-1/m)-1)^(1/n-1)] * Se^(-1/m-1)
    //  Units: [Pa] (since dSw is dimensionless)
    volScalarField computedPcdSw(const volScalarField& Se, const dimensionedScalar Sr) const
    {
        return -(pc0_ / n_) * (1.0 / (1.0 - Sr)) * pow(Se, -1.0 / (n_ - 1.0));
    }
    
    //- Compute relative permeability from effective saturation  
    //  Kr = Se^l * [1 - (1 - Se^(1/m))^m]^2
    //  This is the Mualem-van Genuchten model
    //  Units: [-] (dimensionless)
    volScalarField computeKr(const volScalarField& Se) const
    {
        return pow(Se, m_);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
