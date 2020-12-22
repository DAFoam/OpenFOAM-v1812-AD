/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Test reverse mode matrix vector product dRdWT * Psi

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
    {
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
        tape.setActive();

        // compute dRdWT * psi using reverse mode AD. Here psi is a random vector
        // psi = cos(0.1*cellI) or cos(0.1*(cellI + comp))

        // register input
        forAll(U, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                tape.registerInput(U[cellI][comp]);
            }
        }

        forAll(p, cellI)
        {
            tape.registerInput(p[cellI]);
        }

        forAll(phi, faceI)
        {
            tape.registerInput(phi[faceI]);
        }

        forAll(phi.boundaryField(), patchI)
        {
            forAll(phi.boundaryField()[patchI], faceI)
            {
                tape.registerInput(phi.boundaryFieldRef()[patchI][faceI]);
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

        // register output
        forAll(URes, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                tape.registerOutput(URes[cellI][comp]);
            }
        }

        forAll(pRes, cellI)
        {
            tape.registerOutput(pRes[cellI]);
        }

        forAll(phiRes, faceI)
        {
            tape.registerOutput(phiRes[faceI]);
        }

        forAll(phiRes.boundaryField(), patchI)
        {
            forAll(phiRes.boundaryField()[patchI], faceI)
            {
                tape.registerOutput(phiRes.boundaryFieldRef()[patchI][faceI]);
            }
        }

        tape.setPassive();

        // set seeds
        word productNameSeed = "dRdWTPsi_" + name(myProc) + "_AD_Seeds.txt";
        OFstream fOutSeed(productNameSeed);
        forAll(URes, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar randomSeed = cos(0.1 * (cellI + comp));
                fOutSeed << randomSeed << endl;
                URes[cellI][comp].setGradient(randomSeed.getValue());
            }
        }

        forAll(pRes, cellI)
        {
            scalar randomSeed = cos(0.1 * cellI);
            fOutSeed << randomSeed << endl;
            pRes[cellI].setGradient(randomSeed.getValue());
        }

        forAll(phiRes, faceI)
        {
            scalar randomSeed = cos(0.1 * faceI);
            fOutSeed << randomSeed << endl;
            phiRes[faceI].setGradient(randomSeed.getValue());
        }

        forAll(phiRes.boundaryField(), patchI)
        {
            forAll(phiRes.boundaryField()[patchI], faceI)
            {
                scalar randomSeed = cos(0.1 * (patchI + faceI));
                fOutSeed << randomSeed << endl;
                phiRes.boundaryFieldRef()[patchI][faceI].setGradient(randomSeed.getValue());
            }
        }

        tape.evaluate();

        // output the matrix-vector product to files
        word productName = "dRdWTPsi_" + name(myProc) + "_AD_Values.txt";
        OFstream fOut(productName);
        forAll(U, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar val = U[cellI][comp].getGradient();
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
            scalar val = p[cellI].getGradient();
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
            scalar val = phi[faceI].getGradient();
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
                scalar val = phi.boundaryFieldRef()[patchI][faceI].getGradient();
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
