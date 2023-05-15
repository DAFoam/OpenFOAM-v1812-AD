/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Description
    A brute-force reverse AD solver for DASimpleFoam. 
    Objective function: drag
    Design variable: Volume coordinates
    NOTE: this approach uses a lot of memory!!! Don't use more than 1K mesh cells
    with more than 100 steps.    

\*---------------------------------------------------------------------------*/
#include <codi.hpp>
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addOption
    (
        "patchNames",
        "'(wall)'",
        "List of patch names to compute drag"
    );

    argList::addOption
    (
        "dragDir",
        "'(1 0 0)'",
        "Drag direction"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // read options
    List<wordRe> patchNames;
    if (args.optionFound("patchNames"))
    {
        patchNames = wordReList(args.optionLookup("patchNames")());
    }
    else
    {
        Info<<"drag patchNames not set! Exit."<<endl;
        Info<<"Example: DASimpleFoamReverseAD -patchNames '(wall)' "<<endl;
        return 1;
    }
 
    vector dragDir = {1.0, 0.0, 0.0};
    if (args.optionFound("dragDir"))
    {
        scalarList tmpList=args.optionLookup("dragDir")();
        forAll(tmpList,idxI)
        {
            dragDir[idxI] = tmpList[idxI];
        }
    }
    else
    {
        Info<<"Drag not set! Using default (1 0 0)"<<endl;
    }

    // setup AD inputs
    pointField meshPoints = mesh.points();
    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
    tape.setActive();
    forAll(meshPoints, i)
    {
        for (label j = 0; j < 3; j++)
        {
            tape.registerInput(meshPoints[i][j]);
        }
    }
    mesh.movePoints(meshPoints);

    // run simpleFoam
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    // compute drag
    const surfaceVectorField::Boundary& Sfb = mesh.Sf().boundaryField();
    tmp<volSymmTensorField> tdevRhoReff = turbulence->devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReff().boundaryField();
    vector forces = vector::zero;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == "wall")
        {
            // normal force
            vectorField fN = Sfb[patchI]*p.boundaryField()[patchI];
            // tangential force
            vectorField fT = Sfb[patchI] & devRhoReffb[patchI];
            forAll(fT, faceI) forces += fN[faceI] + fT[faceI]; 
        }
    }
    // project drag to the dragDir
    scalar drag = forces & dragDir;
    Info<<"Drag: "<<drag<<endl;
    // register f output
    tape.registerOutput(drag);
    tape.setPassive();
    drag.setGradient(1.0);
    tape.evaluate();

    // save dFdXv to files
    label nProcs = Pstream::nProcs();
    label myProc = Pstream::myProcNo();
    std::ostringstream np("");
    std::ostringstream mp("");
    np<<nProcs;
    mp<<myProc;
    std::string fName="dDragdXv_"+mp.str()+ "_" + np.str() + ".txt";
    OFstream fOut(fName);
    forAll(meshPoints, i)
    {
        for (label j = 0; j < 3; j++)
        {
            fOut<<meshPoints[i][j].getGradient()<<endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
