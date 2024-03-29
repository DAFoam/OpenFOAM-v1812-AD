/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Calculating wall distance field" << endl;
volScalarField y
(
    IOobject
    (
        "y",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, Zero),
    zeroGradientFvPatchScalarField::typeName
);
y.primitiveFieldRef() = wallDist::New(mesh).y().primitiveField();
y.correctBoundaryConditions();


// Set the mean boundary-layer thickness
dimensionedScalar ybl("ybl", dimLength, Zero);

if (args.found("ybl"))
{
    // If the boundary-layer thickness is provided use it
    ybl.value() = args.opt<scalar>("ybl");
}
else if (args.found("Cbl"))
{
    // Calculate boundary layer thickness as Cbl*mean distance to wall
    ybl.value() = gAverage(y)*args.opt<scalar>("Cbl");
}

Info<< "\nCreating boundary-layer for U of thickness "
    << ybl.value() << " m" << nl << endl;

Info<< "Creating mask field" << endl;
volScalarField mask
(
    IOobject
    (
        "mask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero),
    zeroGradientFvPatchScalarField::typeName
);


// ************************************************************************* //
