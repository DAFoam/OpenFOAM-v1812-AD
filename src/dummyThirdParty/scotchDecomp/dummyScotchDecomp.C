/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"You are trying to use scotch but do not have the scotchDecomp library loaded."
"\nThis message is from the dummy scotchDecomp stub library instead.\n"
"\n"
"Please install scotch and make sure that libscotch.so is in your "
"LD_LIBRARY_PATH.\n"
"The scotchDecomp library can then be built in "
"src/parallel/decompose/decompositionMethods/scotchDecomp\n";


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scotchDecomp::graphPath(const polyMesh& unused) const
{}


void Foam::scotchDecomp::check(const int retVal, const char* str)
{}


Foam::label Foam::scotchDecomp::decomposeSerial
(
    const labelUList& adjncy,
    const labelUList& xadj,
    const UList<scalar>& cWeights,
    List<label>& decomp
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp
(
    const dictionary& decompDict
)
:
    metisLikeDecomp("scotch", decompDict)
{}


Foam::scotchDecomp::scotchDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp("scotch", decompDict, regionName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


// ************************************************************************* //
