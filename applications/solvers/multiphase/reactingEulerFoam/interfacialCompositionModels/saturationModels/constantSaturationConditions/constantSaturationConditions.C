/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "constantSaturationConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(constantSaturationConditions, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        constantSaturationConditions,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::constantSaturationConditions::
constantSaturationConditions(const dictionary& dict)
:
    saturationModel(),
    pSat_("pSat", dimPressure, dict),
    Tsat_("Tsat", dimTemperature, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::constantSaturationConditions::
~constantSaturationConditions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantSaturationConditions::pSat
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "pSat",
            T.mesh().time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        T.mesh(),
        pSat_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantSaturationConditions::pSatPrime
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "pSatPrime",
            T.mesh().time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        T.mesh(),
        dimensionedScalar(dimPressure/dimTemperature, Zero)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantSaturationConditions::lnPSat
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "lnPSat",
            T.mesh().time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        T.mesh(),
        dimensionedScalar("lnPSat", dimless, log(pSat_.value()))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantSaturationConditions::Tsat
(
    const volScalarField& p
) const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "Tsat",
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        p.mesh(),
        Tsat_
    );
}


// ************************************************************************* //
