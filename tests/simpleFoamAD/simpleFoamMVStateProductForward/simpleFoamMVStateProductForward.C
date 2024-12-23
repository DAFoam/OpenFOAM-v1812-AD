/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Test forwad mode matrix vector product dRdW * Psi and compare with FD

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

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createFields.H"
#include "initContinuityErrs.H"

    label myProc = Pstream::myProcNo();

    scalar URef = 1.0;
    scalar pRef = 1.0;
    scalar phiRef = 1.0e-3;

    // AD
    {
        // compute dRdW * psi using forward mode AD. Here psi is a random vector
        // psi = sin(0.1*cellI) or sin(0.1*(cellI + comp))

        // set seeds
        word productNameSeed = "dRdWPsi_" + name(myProc) + "_AD_Seeds.txt";
        OFstream fOutSeed(productNameSeed);
        forAll(U, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar randomSeed = URef * sin(0.1 * (cellI + comp));
                fOutSeed << randomSeed << endl;
                U[cellI][comp].setGradient(randomSeed.getValue());
            }
        }

        forAll(p, cellI)
        {
            scalar randomSeed = pRef * sin(0.1 * cellI);
            fOutSeed << randomSeed << endl;
            p[cellI].setGradient(randomSeed.getValue());
        }

        forAll(phi, faceI)
        {
            scalar randomSeed = phiRef * sin(0.1 * faceI);
            fOutSeed << randomSeed << endl;
            phi[faceI].setGradient(randomSeed.getValue());
        }

        forAll(phi.boundaryField(), patchI)
        {
            forAll(phi.boundaryField()[patchI], faceI)
            {
                scalar randomSeed = phiRef * sin(0.1 * (patchI + faceI));
                fOutSeed << randomSeed << endl;
                phi.boundaryFieldRef()[patchI][faceI].setGradient(randomSeed.getValue());
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
        word productName = "dRdWPsi_" + name(myProc) + "_AD_Values.txt";
        OFstream fOut(productName);
        forAll(URes, cellI)
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

        forAll(pRes, cellI)
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

        forAll(phiRes, faceI)
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

        forAll(phiRes.boundaryField(), patchI)
        {
            forAll(phiRes.boundaryField()[patchI], faceI)
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
