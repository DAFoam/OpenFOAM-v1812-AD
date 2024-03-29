/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd
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

#include "radiativeIntensityRay.H"
#include "fvm.H"
#include "fvDOM.H"
#include "constants.H"

using namespace Foam::constant;

const Foam::word
Foam::radiation::radiativeIntensityRay::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::radiativeIntensityRay
(
    const fvDOM& dom,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const absorptionEmissionModel& absorptionEmission,
    const blackBodyEmission& blackBody,
    const label rayId
)
:
    dom_(dom),
    mesh_(mesh),
    absorptionEmission_(absorptionEmission),
    blackBody_(blackBody),
    I_
    (
        IOobject
        (
            "I" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qr_
    (
        IOobject
        (
            "qr" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qin_
    (
        IOobject
        (
            "qin" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qem_
    (
        IOobject
        (
            "qem" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    d_(Zero),
    dAve_(Zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    myRayId_(rayId)
{
    scalar sinTheta = sin(theta);
    scalar cosTheta = cos(theta);
    scalar sinPhi = sin(phi);
    scalar cosPhi = cos(phi);

    omega_ = 2.0*sinTheta*sin(deltaTheta/2.0)*deltaPhi;
    d_ = vector(sinTheta*sinPhi, sinTheta*cosPhi, cosTheta);
    dAve_ = vector
    (
        sinPhi
       *sin(0.5*deltaPhi)
       *(deltaTheta - cos(2.0*theta)
       *sin(deltaTheta)),
        cosPhi
       *sin(0.5*deltaPhi)
       *(deltaTheta - cos(2.0*theta)
       *sin(deltaTheta)),
        0.5*deltaPhi*sin(2.0*theta)*sin(deltaTheta)
    );

    if (mesh_.nSolutionD() == 2)
    {
        // Omega for 2D
        omega_ = deltaPhi;

        // dAve for 2D
        dAve_ = vector
        (
            2*sinPhi*sin(0.5*deltaPhi),
            2*cosPhi*sin(0.5*deltaPhi),
            0
        );

        vector meshDir(Zero);
        if (dom_.meshOrientation() != vector::zero)
        {
            meshDir = dom_.meshOrientation();
        }
        else
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh_.geometricD()[cmpt] == -1)
                {
                    meshDir[cmpt] = 1;
                }
            }
        }
        const vector normal(vector(0, 0, 1));

        const tensor coordRot = rotationTensor(normal, meshDir);

        dAve_ = coordRot & dAve_;
        d_ = coordRot & d_;

    }
    else if (mesh_.nSolutionD() == 1)
    {
        vector meshDir(Zero);
        if (dom_.meshOrientation() != vector::zero)
        {
            meshDir = dom_.meshOrientation();
        }
        else
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh_.geometricD()[cmpt] == 1)
                {
                    meshDir[cmpt] = 1;
                }
            }
        }
        const vector normal(vector(1, 0, 0));

        dAve_ = (dAve_ & normal)*meshDir;
        d_ = (d_ & normal)*meshDir;

        // Omega normalization for 1D
        omega_ /= 2;
    }

    autoPtr<volScalarField> IDefaultPtr;

    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(rayId) + "_" + name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // Check if field exists and can be read
        if (IHeader.typeHeaderOk<volScalarField>(true))
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IDefault field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IDefault",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }

            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::~radiativeIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::radiativeIntensityRay::correct()
{
    // Reset boundary heat flux to zero
    qr_.boundaryFieldRef() = 0.0;

    scalar maxResidual = -GREAT;

    forAll(ILambda_, lambdaI)
    {
        const volScalarField& k = dom_.aLambda(lambdaI);

        const surfaceScalarField Ji(dAve_ & mesh_.Sf());

        fvScalarMatrix IiEq
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
          + fvm::Sp(k*omega_, ILambda_[lambdaI])
        ==
            1.0/constant::mathematical::pi*omega_
           *(
                (k - absorptionEmission_.aDisp(lambdaI))
               *blackBody_.bLambda(lambdaI)

              + absorptionEmission_.E(lambdaI)/4
            )
        );

        IiEq.relax();

        const solverPerformance ILambdaSol = solve
        (
            IiEq,
            mesh_.solver("Ii")
        );

        const scalar initialRes =
            ILambdaSol.initialResidual()*omega_/dom_.omegaMax();

        maxResidual = max(initialRes, maxResidual);
    }

    return maxResidual;
}


void Foam::radiation::radiativeIntensityRay::addIntensity()
{
    I_ = dimensionedScalar(dimMass/pow3(dimTime), Zero);

    forAll(ILambda_, lambdaI)
    {
        I_ += ILambda_[lambdaI];
    }
}


// ************************************************************************* //
