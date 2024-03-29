/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fv.H"
#include "objectRegistry.H"
#include "solution.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fv::gradScheme<Type>> Foam::fv::gradScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing gradScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Grad scheme not specified" << endl << endl
            << "Valid grad schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto cstrIter = IstreamConstructorTablePtr_->cfind(schemeName);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(schemeData)
            << "Unknown grad scheme "
            << schemeName << nl << nl
            << "Valid grad schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gradScheme<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    if (!this->mesh().changing() && this->mesh().cache(name))
    {
        if (!mesh().objectRegistry::template foundObject<GradFieldType>(name))
        {
            solution::cachePrintMessage("Calculating and caching", name, vsf);
            tmp<GradFieldType> tgGrad = calcGrad(vsf, name);
            regIOobject::store(tgGrad.ptr());
        }

        solution::cachePrintMessage("Retrieving", name, vsf);
        GradFieldType& gGrad =
            mesh().objectRegistry::template lookupObjectRef<GradFieldType>
            (
                name
            );

        if (gGrad.upToDate(vsf))
        {
            return gGrad;
        }
        else
        {
            solution::cachePrintMessage("Deleting", name, vsf);
            gGrad.release();
            delete &gGrad;

            solution::cachePrintMessage("Recalculating", name, vsf);
            tmp<GradFieldType> tgGrad = calcGrad(vsf, name);

            solution::cachePrintMessage("Storing", name, vsf);
            regIOobject::store(tgGrad.ptr());
            GradFieldType& gGrad =
                mesh().objectRegistry::template lookupObjectRef<GradFieldType>
                (
                    name
                );

            return gGrad;
        }
    }
    else
    {
        if (mesh().objectRegistry::template foundObject<GradFieldType>(name))
        {
            GradFieldType& gGrad =
                mesh().objectRegistry::template lookupObjectRef<GradFieldType>
                (
                    name
                );

            if (gGrad.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vsf);
                gGrad.release();
                delete &gGrad;
            }
        }

        solution::cachePrintMessage("Calculating", name, vsf);
        return calcGrad(vsf, name);
    }
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gradScheme<Type>::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf
) const
{
    return grad(vsf, "grad(" + vsf.name() + ')');
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::gradScheme<Type>::grad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    tmp<GradFieldType> tgrad = grad(tvsf());
    tvsf.clear();
    return tgrad;
}


// ************************************************************************* //
