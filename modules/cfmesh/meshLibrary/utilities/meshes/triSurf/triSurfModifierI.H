/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::pointField& Foam::Module::triSurfModifier::pointsAccess()
{
    return const_cast<pointField&>(surface_.points());
}


inline Foam::Module::LongList<Foam::labelledTri>&
Foam::Module::triSurfModifier::facetsAccess()
{
    return const_cast<LongList<labelledTri>&>(surface_.facets());
}


inline Foam::Module::edgeLongList&
Foam::Module::triSurfModifier::featureEdgesAccess()
{
    return const_cast<edgeLongList&>(surface_.featureEdges());
}


inline Foam::geometricSurfacePatchList&
Foam::Module::triSurfModifier::patchesAccess()
{
    return const_cast<geometricSurfacePatchList&>(surface_.patches());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
