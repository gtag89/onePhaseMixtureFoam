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

References:
    van Genuchten, M.Th. (1980). A closed-form equation for predicting the 
    hydraulic conductivity of unsaturated soils. Soil Science Society of 
    America Journal, 44(5), 892-898.
*/

class vanGenuchten
{
private:
    
    // Van Genuchten shape parameter m [-]
    dimensionedScalar m_;
    
    // Van Genuchten shape parameter n [-]
    dimensionedScalar n_;
    
    // Pore connectivity parameter l [-]
    dimensionedScalar l_;
    
    // Air-entry pressure field pc0 [Pa]
    volScalarField pc0_;

public:

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
    
    // Compute capillary pressure from effective saturation
    volScalarField computePc(const volScalarField& Se) const
    {
        return pc0_ * pow(pow(Se, -1.0/m_) - 1.0, 1.0/n_);
    }
    
    // Compute derivative of capillary pressure with respect to water saturation
    volScalarField computedPcdSw(const volScalarField& Se, const dimensionedScalar Sr) const
    {
        return -(pc0_ / (m_ * n_)) * (1.0 / (1.0 - Sr)) * pow(pow(Se, -1.0/m_) 
               - 1.0, 1.0/n_ - 1.0) * pow(Se, -1.0/m_ - 1.0);
    }
    
    // Compute relative permeability from effective saturation  
    volScalarField computeKr(const volScalarField& Se) const
    {
        return pow(Se, 0.5) * pow(1 - pow(1 - pow(Se, 1.0/m_), m_), 2);
    }
};

#endif

// ************************************************************************* //
