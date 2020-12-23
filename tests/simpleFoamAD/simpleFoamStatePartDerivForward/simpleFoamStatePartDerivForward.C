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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{

    argList::addOption
    (
        "cell",
        "152",
        "Which cell to compute dRdW"
    );

    argList::addOption
    (
        "proc",
        "0",
        "which procesor to compute"
    );

    argList::addOption
    (
        "var",
        "U",
        "Which residual to compute"
    );

    argList::addOption
    (
        "comp",
        "0",
        "Which component of residual to compute (for U only)"
    );

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createFields.H"

    label cellI = readLabel(args.optionLookup("cell")());
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
    word var = word(args.optionLookup("var")());

    label myProc = Pstream::myProcNo();
    {
        // compute a row of dRdW using forward mode AD.

        // register input
        if (procI == myProc)
        {
            if ("U" == var)
            {
                U[cellI][comp].setGradient(1.0);
            }
            else if ("p" == var)
            {
                p[cellI].setGradient(1.0);
            }
            else if ("phi" == var)
            {
                phi[cellI].setGradient(1.0);
            }

        }

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // URes
        fvVectorMatrix UEqn(fvm::div(phi, U) + turbulence->divDevReff(U));
        UEqn.relax();
        volVectorField URes = (UEqn & U) + fvc::grad(p);
        // pRes
        label pRefCell = 0;
        scalar pRefValue = 0.0;
        volScalarField rAU(1.0 / UEqn.A());
        volVectorField HbyA("HbyA", U);
        HbyA = rAU * UEqn.H();
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, p);
        fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));
        pEqn.setReference(pRefCell, pRefValue);
        volScalarField pRes = pEqn & p;
        // phiRes
        surfaceScalarField phiRes = phiHbyA - pEqn.flux() - phi;
        
        // output the matrix-vector product to files
        word productName = "dRdW_" + var + name(comp) + "_Idx_" + name(cellI) + "_ProcI"  + name(myProc) + "_FAD.txt";
        OFstream fOut(productName);
        forAll(U, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar val = URes[cellI][comp].getGradient();
                if (fabs(val) > 1e-16)
                {
                    fOut << val << endl;
                }
                else
                {
                    fOut << "0" << endl;
                }
            }
        }

        forAll(p, cellI)
        {
            scalar val = pRes[cellI].getGradient();
            if (fabs(val) > 1e-16)
            {
                fOut << val << endl;
            }
            else
            {
                fOut << "0" << endl;
            }
        }

        forAll(phi, faceI)
        {
            scalar val = phiRes[faceI].getGradient();
            if (fabs(val) > 1e-16)
            {
                fOut << val << endl;
            }
            else
            {
                fOut << "0" << endl;
            }
        }

        forAll(phi.boundaryField(), patchI)
        {
            forAll(phi.boundaryField()[patchI], faceI)
            {
                scalar val = phiRes.boundaryFieldRef()[patchI][faceI].getGradient();
                if (fabs(val) > 1e-16)
                {
                    fOut << val << endl;
                }
                else
                {
                    fOut << "0" << endl;
                }
            }
        }
    }

    Info << "Done!" << endl;
    return 0;
}

// ************************************************************************* //
