/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "fieldCoordinateSystemTransform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldCoordinateSystemTransform, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fieldCoordinateSystemTransform,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
fieldCoordinateSystemTransform
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_),
    csysPtr_
    (
        coordinateSystem::New(mesh_, dict, coordinateSystem::typeName_())
    )
{
    read(dict);

    Info<< type() << " " << name << ":" << nl
        << "   Applying " << (csysPtr_->uniform() ? "" : "non-")
        << "uniform transformation from global Cartesian to local "
        << *csysPtr_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word
Foam::functionObjects::fieldCoordinateSystemTransform::transformFieldName
(
    const word& fieldName
) const
{
    return fieldName + ":Transformed";
}


const Foam::surfaceTensorField&
Foam::functionObjects::fieldCoordinateSystemTransform::srotTensor() const
{
    typedef surfaceTensorField FieldType;
    typedef surfaceTensorField::Boundary BoundaryType;

    if (!rotTensorSurface_.valid())
    {
        tensorField rotations(csysPtr_->R(mesh_.faceCentres()));

        rotTensorSurface_.reset
        (
            new FieldType
            (
                IOobject
                (
                    "surfRotation",
                    mesh_.objectRegistry::instance(),
                    mesh_.objectRegistry::db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false // no register
                ),
                mesh_,
                dimless,
                std::move(rotations)
                // calculatedType
            )
        );

        auto& rot = *rotTensorSurface_;

        // Boundaries
        BoundaryType& bf = const_cast<BoundaryType&>(rot.boundaryField());

        forAll(bf, patchi)
        {
            bf[patchi] = csysPtr_->R(bf[patchi].patch().patch().faceCentres());
        }
    }

    return *rotTensorSurface_;
}


const Foam::volTensorField&
Foam::functionObjects::fieldCoordinateSystemTransform::vrotTensor() const
{
    typedef volTensorField FieldType;
    typedef volTensorField::Boundary BoundaryType;

    if (!rotTensorVolume_.valid())
    {
        tensorField rotations(csysPtr_->R(mesh_.cellCentres()));

        rotTensorVolume_.reset
        (
            new FieldType
            (
                IOobject
                (
                    "volRotation",
                    mesh_.objectRegistry::instance(),
                    mesh_.objectRegistry::db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false // no register
                ),
                mesh_,
                dimless,
                std::move(rotations)
                // calculatedType
            )
        );

        auto& rot = *rotTensorVolume_;

        // Boundaries
        BoundaryType& bf = const_cast<BoundaryType&>(rot.boundaryField());

        forAll(bf, patchi)
        {
            bf[patchi] = csysPtr_->R(bf[patchi].patch().patch().faceCentres());
        }
    }

    return *rotTensorVolume_;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::read
(
    const dictionary& dict
)
{
    if (fvMeshFunctionObject::read(dict))
    {
        fieldSet_.read(dict);
        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::execute()
{
    fieldSet_.updateSelection();

    for (const word& fieldName : fieldSet_.selection())
    {
        transform<scalar>(fieldName);
        transform<vector>(fieldName);
        transform<sphericalTensor>(fieldName);
        transform<symmTensor>(fieldName);
        transform<tensor>(fieldName);
    }

    // Finished with these
    rotTensorSurface_.clear();
    rotTensorVolume_.clear();

    return true;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::write()
{
    forAllConstIters(fieldSet_, iter)
    {
        writeObject(transformFieldName(iter()));
    }

    return true;
}


// ************************************************************************* //
