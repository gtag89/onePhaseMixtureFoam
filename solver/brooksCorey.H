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

References:
    Brooks, R.H. and Corey, A.T. (1964). Hydraulic properties of porous media.
    Hydrology Papers, Colorado State University, Fort Collins, Colorado.
*/

class brooksCorey
{
private:

    // Pore size distribution index lambda [-]
    dimensionedScalar lambda_;
    
    // Relative permeability exponent eta [-]
    dimensionedScalar eta_;
    
    // Displacement (air-entry) pressure field pd [Pa]
    volScalarField pcd_;

public:

    brooksCorey
    (
        const dimensionedScalar lambda,
        const dimensionedScalar eta, 
        const volScalarField& pcd
    )
    :
        lambda_(lambda), 
        eta_(eta), 
        pcd_(pcd)
    {}

    // Compute capillary pressure from effective saturation
    volScalarField computePc(const volScalarField& Se) const
    {
        return pcd_ * pow(Se, -1.0/lambda_);
    }
    
    // Compute derivative of capillary pressure with respect to water saturation
    volScalarField computedPcdSw(const volScalarField& Se, const dimensionedScalar Sr) const
    {
        return -(pcd_ / lambda_) * (1.0 / (1.0 - Sr)) * pow(Se, -1.0 / (lambda_ - 1.0));
    }
    
    // Compute relative permeability from effective saturation  
    volScalarField computeKr(const volScalarField& Se) const
    {
        return pow(Se, eta_);
    }
};

#endif

// ************************************************************************* //
