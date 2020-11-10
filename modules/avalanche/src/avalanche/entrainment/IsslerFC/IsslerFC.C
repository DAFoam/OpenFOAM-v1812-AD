/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | faSavageHutterFOAM
    \\  /    A nd           | Copyright (C) 2017 Matthias Rauter
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

Author
    Matthias Rauter matthias.rauter@uibk.ac.at

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "IsslerFC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentModels
{
    defineTypeNameAndDebug(IsslerFC, 0);
    addToRunTimeSelectionTable(entrainmentModel, IsslerFC, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentModels::IsslerFC::IsslerFC
(
    const dictionary& entrainmentProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& hentrain,
    const areaScalarField& pb,
    const areaVectorField& tau
)
:
    entrainmentModel(type(), entrainmentProperties, Us, h, hentrain, pb, tau),
    tauc_("tauc", coeffDict_),
    mu_("mu", coeffDict_),
    K_("K", coeffDict_),
    gs_(Us.db().lookupObject<areaVectorField>("gs")),
    gn_(Us.db().lookupObject<areaScalarField>("gn"))
{
    Info << "    " << tauc_ << endl << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::entrainmentModels::IsslerFC::Sm() const
{
    const dimensionedScalar smallVel("small", dimVelocity, SMALL);

    const areaScalarField u(mag(Us_));
    const areaScalarField gs(mag(gs_));
    const areaScalarField gamma_c(5./2.*u/h_);

    const areaScalarField uthr
    (
        h_*gamma_c*gamma_c/2.
      * K_*gamma_c
      / ( h_*(gs+mu_*gn_)+4.*K_*gamma_c*gamma_c )
    );

    const areaScalarField weinf
    (
        (h_*(gs+mu_*gn_)+4*K_*pow(gamma_c,2))
      / (h_*gamma_c+smallVel)
    );

    Sm_ = pos(u-uthr)*weinf*(1-uthr/(u+smallVel));

    const areaScalarField hlimit(h_*mag(gs_)-(tauc_-5*K_*gamma_c*gamma_c));
    const areaScalarField hlimit2(h_-dimensionedScalar("small", dimLength, 1e-2));

    Sm_  = pos(hlimit2)*pos(hlimit)*Sm_;

    Sm_ = min(Sm_, hentrain_/Us_.db().time().deltaT());

    Info << "IsslerFC:min/max(Sm) =  " << min(Sm_) << " / " << max(Sm_) << endl;

    return Sm_;
}


bool Foam::entrainmentModels::IsslerFC::read
(
    const dictionary& entrainmentProperties
)
{
    readDict(type(), entrainmentProperties);

    coeffDict_.readEntry("tauc", tauc_);

    return true;
}


// ************************************************************************* //
