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

    Verify brute force automatic differentiation against finite-difference.
    The design variable is the velocity at the inlet and the objective function
    is the sum of x velocity of all cells.  

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar calcObj(const volVectorField& U)
{
    scalar obj = scalar(0.0);
    forAll(U, idxI)
    {
        obj += U[idxI][0];
    }

    return obj;
}

void setUIn(const fvMesh& mesh, const word& patchName, volVectorField& U, const scalar& UBC)
{
    label patchI = mesh.boundaryMesh().findPatchID(patchName);

    forAll(U.boundaryField()[patchI], faceI)
    {
        U.boundaryFieldRef()[patchI][faceI][0] = UBC;
    }
}

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "mode",
        "AD",
        "Options are: AD (foward-mode or reverse-mode AD) or FD (finite-difference)" 
    );

    argList::addOption
    (
        "UIn",
        "10.0",
        "Velocity at inlet"
    );

    //#include "postProcess.H"

    //#include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    word mode="AD";
    if (args.optionFound("mode")) 
    {
        mode = word(args.optionLookup("mode")());
    }
    else 
    {
        Info<<"mode arg not found! Use the default -mode AD"<<endl
            <<"Options are:"<<endl
            <<"AD (foward-mode or reverse-mode AD)"<<endl
            <<"FD (finite-difference)"<<endl 
            <<"Example: simpleROMFoam -mode AD"<<endl;
    }

    scalar UxIn = 10.0;
    if (args.optionFound("UIn"))
    {
        UxIn = readScalar(args.optionLookup("UIn")());
    }
    else
    {
        Info<<"UIn arg not found! Use the default -UIn 10.0"<<endl;
    }

    if (mode == "AD")
    {
#ifdef CODI_AD_FORWARD
        UxIn.setGradient(1.0);
#endif
#ifdef CODI_AD_REVERSE
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
        tape.setActive();
        tape.registerInput(UxIn);
#endif
    }

    setUIn(mesh, "inlet", U, UxIn);    

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
        Info<<"Obj: "<<calcObj(U)<<endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    scalar objFunc = calcObj(U);

    if (mode == "AD")
    {
#ifdef CODI_AD_FORWARD
        Info<<"f: "<<objFunc<<endl;
        Info<<"df/dx: "<< objFunc.getGradient()<<endl;
#endif
#ifdef CODI_AD_REVERSE
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
        tape.registerOutput(objFunc);
        tape.setPassive();
        objFunc.setGradient(1.0);
        tape.evaluate();
        Info<<"f: "<<objFunc<<endl;
        Info<<"df/dx: "<< UxIn.getGradient()<<endl;
#endif
    }
    else
    {
        Info<<"f: "<<objFunc<<endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
