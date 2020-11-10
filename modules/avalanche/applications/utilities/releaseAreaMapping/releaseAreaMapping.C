/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | faSavageHutterFOAM
    \\  /    A nd           | Copyright (C) 2017 Matthias Rauter
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    releaseAreaMapping

Description
    Initialization of finite area fields,
    used to create release areas for avalanche simulations.

Author
    Matthias Rauter matthias.rauter@uibk.ac.at

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "HormannAgathos.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "(avalanche)\n"
        "Initialization of finite area fields,"
        " used to create release areas for avalanche simulations."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFaMesh.H"

    areaVectorField c = aMesh.areaCentres();
    areaVectorField n = aMesh.faceAreaNormals();


    IOdictionary releaseArea
    (
         IOobject
         (
              "releaseArea",
              runTime.constant(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE
         )
    );

    const PtrList<entry> fields
    (
        releaseArea.lookup("fields")
    );

    wordList fieldNames;
    fieldNames.setSize(fields.size());

    forAll(fieldNames, fieldsI)
    {
        const entry& fieldsInfo = fields[fieldsI];

        if (!fieldsInfo.isDict())
        {
            FatalIOErrorIn("releaseAreaMapping.C", releaseArea)
                << "Entry " << fieldsInfo << " in fields section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }
        fieldNames[fieldsI] = fieldsInfo.keyword();

        dictionary fieldsDict = fieldsInfo.dict();

        Info<< "Reading field " << fieldNames[fieldsI] << endl;

        areaScalarField f
        (
            IOobject
            (
                fieldNames[fieldsI],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh
        );

        dimensionedScalar fdefault = fieldsDict.lookupOrDefault<dimensionedScalar>("default", dimensionedScalar("default", f.dimensions(), -1e100));

        if (fdefault.value() > -1e99)
        {
            Info<< "Setting field " << fieldNames[fieldsI]
                << " to default value" << endl;
            f.ref() = fdefault;
        }
        string str = fieldNames[fieldsI];
        str[0] = toupper(str[0]);

        Info<< "Creating field " << str << endl;

        volScalarField F
        (
            IOobject
            (
                str,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("0", f.dimensions(), 0)
        );

        // Create volume-to surface mapping object
        volSurfaceMapping vsm(aMesh);


        Info<< "Reading regions" << endl;

        const PtrList<entry> regions
        (
            fieldsDict.lookup("regions")
        );

        wordList areaNames;
        areaNames.setSize(regions.size());

        forAll(regions, areaI)
        {
            const entry& regionInfo = regions[areaI];

            Info<< "processing region " << regions[areaI].keyword() << endl;
            if (!regionInfo.isDict())
            {
                FatalIOErrorIn("releaseAreaMapping.C", releaseArea)
                    << "Entry " << regionInfo << " in boundary section is not a"
                    << " valid dictionary." << exit(FatalIOError);
            }
            areaNames[areaI] = regionInfo.keyword();

            dictionary areaDict = regionInfo.dict();

            word type;
            vector offset;
            List<point2D> points;

            areaDict.lookup("type") >> type;


            if (type == "polygon")
            {
                scalar finit;
                List<vector> vertices;

                areaDict.lookup("offset") >> offset;
                areaDict.lookup("vertices") >> vertices;
                areaDict.lookup("value") >> finit;
                Switch projectNormal = areaDict.lookupOrDefault<Switch>("projectToNormal", false);

                areaScalarField projection
                (
                    IOobject
                    (
                        "projection",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    aMesh,
                    dimensionedScalar("0", dimless, 1)
                );

                if (projectNormal)
                {
                    projection = n&vector(0, 0, -1);
                }

                points.resize(vertices.size());

                forAll(vertices, vI)
                {
                    points[vI] = point2D(vertices[vI].x()+offset.x(), vertices[vI].y()+offset.y());
                }
                HormannAgathos polygon(points, 0.001);

                forAll(c.internalField(), i)
                {
                    if (polygon.evaluate(point2D(c[i].x(), c[i].y())) != HormannAgathos::POINT_OUTSIDE)
                    {
                        f[i] = finit*projection[i];
                    }
                }
            }
            else if (type == "polygonlinear")
            {
                scalar finit;
                scalar fx, fy, fz;
                scalar x0, y0, z0;
                List<vector> vertices;

                areaDict.lookup("offset") >> offset;
                areaDict.lookup("vertices") >> vertices;
                areaDict.lookup("valueAtZero") >> finit;
                x0 = areaDict.lookupOrDefault<scalar>("x0", 0);
                y0 = areaDict.lookupOrDefault<scalar>("y0", 0);
                z0 = areaDict.lookupOrDefault<scalar>("z0", 0);
                fx = areaDict.lookupOrDefault<scalar>("dfdx", 0);
                fy = areaDict.lookupOrDefault<scalar>("dfdy", 0);
                fz = areaDict.lookupOrDefault<scalar>("dfdz", 0);
                Switch projectNormal = areaDict.lookupOrDefault<Switch>("projectToNormal", false);

                areaScalarField projection
                (
                    IOobject
                    (
                        "projection",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    aMesh,
                    dimensionedScalar("0", dimless, 1)
                );

                if (projectNormal)
                {
                    projection = n&vector(0, 0, -1);
                }

                Info<< "linear function with" << nl
                    << "     x0 = " << x0 << ", dfdx = " << fx << nl
                    << "     y0 = " << y0 << ", dfdy = " << fy << nl
                    << "     z0 = " << z0 << ", dfdz = " << fz << endl;

                points.resize(vertices.size());

                forAll(vertices, vI)
                {
                    points[vI] = point2D(vertices[vI].x()+offset.x(), vertices[vI].y()+offset.y());
                }
                HormannAgathos polygon(points, 0.001);

                forAll(c.internalField(), i)
                {
                    if (polygon.evaluate(point2D(c[i].x(), c[i].y())) != HormannAgathos::POINT_OUTSIDE)
                    {
                        f[i] = (finit + fx*(c[i].x()-x0) + fy*(c[i].y()-y0) + fz*(c[i].z()-z0))*projection[i];
                    }
                }
            }
            else if (type == "sphere")
            {
                vector center;
                scalar rad;
                scalar scale;

                areaDict.lookup("center") >> center;
                areaDict.lookup("r") >> rad;
                scale = areaDict.lookupOrDefault<scalar>("scale", scalar(1));

                forAll(c.internalField(), i)
                {
                    vector c_cs = c[i]-center;
                    scalar c_csmag = Foam::mag(c_cs);
                    if (c_csmag < rad)
                    {
                        scalar cy = -n[i] & c_cs;
                        scalar cx = Foam::sqrt(Foam::sqr(Foam::mag(c_cs))-Foam::sqr(cy));
                        f[i] = scale*(Foam::sqrt(Foam::sqr(rad)-Foam::sqr(cx))-cy);
                    }
                }
            }
        }

        Info<< "Writing fields" << endl;

        vsm.mapToVolume(f, F.boundaryFieldRef());

        f.write();
        F.write();
    }

    Info<< nl << "End" << endl;
    return 0;
}


// ************************************************************************* //
