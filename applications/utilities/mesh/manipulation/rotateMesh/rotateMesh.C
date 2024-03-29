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

Application
    rotateMesh

Group
    grpMeshManipulationUtilities

Description
    Rotates the mesh and fields from the direction n1 to direction n2.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "transformGeometricField.H"
#include "IOobjectList.H"

using namespace Foam;

template<class GeometricField>
void RotateFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const tensor& rotT
)
{
    // Objects of field type
    IOobjectList fields(objects.lookupClass(GeometricField::typeName));

    forAllConstIters(fields, fieldIter)
    {
        Info<< "    Rotating " << fieldIter()->name() << endl;

        GeometricField fld(*fieldIter(), mesh);
        transform(fld, dimensionedTensor(rotT), fld);
        fld.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Rotate mesh points and vector/tensor fields\n"
        "Rotation from the <from> vector to the <to> vector"
    );

    timeSelector::addOptions();

    argList::addArgument("from", "The vector to rotate from");
    argList::addArgument("to",   "The vector to rotate to");

    #include "setRootCase.H"
    #include "createTime.H"

    const vector n1(args.get<vector>(1).normalise());
    const vector n2(args.get<vector>(2).normalise());

    const tensor rotT(rotationTensor(n1, n2));

    {
        pointIOField points
        (
            IOobject
            (
                "points",
                runTime.findInstance(polyMesh::meshSubDir, "points"),
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        points = transform(rotT, points);

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        Info<< "Writing points into directory " << points.path() << nl << endl;
        points.write();
    }


    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        RotateFields<volVectorField>(mesh, objects, rotT);
        RotateFields<volSphericalTensorField>(mesh, objects, rotT);
        RotateFields<volSymmTensorField>(mesh, objects, rotT);
        RotateFields<volTensorField>(mesh, objects, rotT);

        RotateFields<surfaceVectorField>(mesh, objects, rotT);
        RotateFields<surfaceSphericalTensorField>(mesh, objects, rotT);
        RotateFields<surfaceSymmTensorField>(mesh, objects, rotT);
        RotateFields<surfaceTensorField>(mesh, objects, rotT);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
