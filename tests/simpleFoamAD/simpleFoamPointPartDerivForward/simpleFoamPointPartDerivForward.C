/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Test forward mode partial derivative dRdW

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "OFstream.H"
#include "IOstreamOption.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{

    argList::addOption
    (
        "point",
        "152",
        "Which point to perturb"
    );

    argList::addOption
    (
        "comp",
        "2",
        "Which component to perturb"
    );

    argList::addOption
    (
        "proc",
        "0",
        "which procesor to compute"
    );

    argList::addOption(
        "mode",
        "FD",
        "Which mode to run");

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createFields.H"

    label pointI = readLabel(args.optionLookup("point")());
    label procI = 0;
    if (args.optionFound("proc"))
    {
        procI = readLabel(args.optionLookup("proc")());
    }
    label comp = 0;
    if (args.optionFound("comp"))
    {
        comp = readLabel(args.optionLookup("comp")());
    }
    word mode = "FD";
    if (args.optionFound("mode"))
    {
        mode = word(args.optionLookup("mode")());
        Info << "Using mode: " << mode << endl;
    }
    else
    {
        Info << "Using default mode: " << mode << endl;
    }


    label myProc = Pstream::myProcNo();
    if(mode=="AD")
    {
        // compute dFdXv using forward mode AD.

        pointField meshPoints = mesh.points();
        // register input
        if (procI == myProc)
        {
            meshPoints[pointI][comp].setGradient(1.0);
        }
        mesh.movePoints(meshPoints);

        volVectorField gradP = fvc::grad(p);
        scalar f = 0.0;

        if(Pstream::master())
        {
            forAll(gradP, idxI) f+=mag(gradP[idxI]);
        }

        word productName = "dFdXv_ProcI"  + name(myProc) + "_FAD.txt";
        std::ofstream fOut;
        fOut.open(productName, std::ios_base::app); // append instead of overwrite
        scalar deriv = f.getGradient();
        if (fabs(deriv) > 1e-16)
        {
            fOut << deriv << std::endl;
        }
        else
        {
            fOut << "0" << std::endl;
        }

    }
    if (mode=="FD")
    {
        scalar eps = 1e-6;
        // compute dFdW using FD.
        volVectorField gradPRef = fvc::grad(p);
        scalar fRef = 0.0;
        
        if(Pstream::master())
        {   
            forAll(gradPRef, idxI) fRef+=mag(gradPRef[idxI]);
        }

        // perturb
        pointField meshPoints = mesh.points();
        // register input
        if (procI == myProc)
        {
            meshPoints[pointI][comp] += eps;
        }
        mesh.movePoints(meshPoints);

        volVectorField gradP = fvc::grad(p);
        scalar f = 0.0;

        if(Pstream::master())
        {
            forAll(gradP, idxI) f+=mag(gradP[idxI]);
        }

        // output the matrix-vector product to files
        word productName = "dFdXv_ProcI"  + name(myProc) + "_FD.txt";
        std::ofstream fOut;
        fOut.open(productName, std::ios_base::app); // append instead of overwrite
        scalar deriv = (f-fRef)/eps;
        //fOut << "pointI: " << pointI << " comp: "<< comp << " val: ";
        if (fabs(deriv) > 1e-16)
        {
            fOut << deriv << std::endl;
        }
        else
        {
            fOut << "0" << std::endl;
        }

    }


    Info << "Done!" << endl;
    return 0;
}

// ************************************************************************* //
