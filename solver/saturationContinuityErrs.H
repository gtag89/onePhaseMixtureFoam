/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

Global
    saturationContinuityError

Description
    Calculates and prints the saturation continuity errors for unsaturated
    flow simulations. This monitors conservation of the saturation field
    to detect numerical issues or mass balance problems.
    
    The errors calculated are:
    - sumLocalContErr: Sum of absolute local saturation changes (always positive)
    - globalContErr: Net global saturation change (can be positive or negative)  
    - cumulativeContErr: Running total of global errors over time

\*---------------------------------------------------------------------------*/

{
    // ==========================================================================
    // SATURATION MASS CONSERVATION CHECK
    // ==========================================================================
    
    // Calculate total saturation in domain at previous time step
    // This serves as the normalization factor for relative errors
    dimensionedScalar totalMass = fvc::domainIntegrate(SwPrev);
    
    // Prevent division by zero if domain is completely dry
    if (totalMass.value() < SMALL)
    {
        totalMass.value() = SMALL;
        WarningInFunction
            << "Very small total saturation detected: " << totalMass.value() << nl
            << "Continuity error calculation may be inaccurate" << endl;
    }
    
    // Sum of absolute local saturation changes (normalized)
    // This represents the total magnitude of saturation redistribution
    scalar sumLocalContErr = 
        (fvc::domainIntegrate(mag(Sw - SwPrev)) / totalMass).value();
    
    // Net global saturation change (normalized)
    // This should be zero for true conservation (no sources/sinks)
    scalar globalContErr = 
        (fvc::domainIntegrate(Sw - SwPrev) / totalMass).value();
    
    // Accumulate global error over simulation time
    cumulativeContErr += globalContErr;
    
    // ==========================================================================
    // OUTPUT CONTINUITY ERROR INFORMATION
    // ==========================================================================
    Info<< "Saturation continuity errors:" << nl
        << "    Sum local     = " << sumLocalContErr << nl
        << "    Global        = " << globalContErr << nl  
        << "    Cumulative    = " << cumulativeContErr << nl;
}

// ************************************************************************* //
