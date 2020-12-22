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
    // FD
    {
        // compute dRdW * psi using finite-difference. Here psi is a random vector
        // psi = sin(0.1*cellI) or sin(0.1*(cellI + comp))

        // FD perturbation
        scalar eps = 1.0e-6;

        // ref Res
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
        // perturb states
        forAll(U, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar randomSeed = URef * Foam::sin(0.1 * (cellI + comp));
                U[cellI][comp] += randomSeed * eps;
            }
        }

        forAll(p, cellI)
        {
            scalar randomSeed = pRef * Foam::sin(0.1 * cellI);
            p[cellI] += randomSeed * eps;
        }

        forAll(phi, faceI)
        {
            scalar randomSeed = phiRef * Foam::sin(0.1 * faceI);
            phi[faceI] += randomSeed * eps;
        }

        forAll(phi.boundaryField(), patchI)
        {
            forAll(phi.boundaryField()[patchI], faceI)
            {
                scalar randomSeed = phiRef * Foam::sin(0.1 * (patchI + faceI));
                phi.boundaryFieldRef()[patchI][faceI] += randomSeed * eps;
            }
        }

        // compute perturbed Res
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();
        fvVectorMatrix UEqnP(fvm::div(phi, U) + turbulence->divDevReff(U));
        UEqnP.relax();
        volVectorField UResP = (UEqnP & U) + fvc::grad(p);
        volScalarField rAUP(1.0 / UEqnP.A());
        volVectorField HbyAP("HbyA", U);
        HbyAP = rAUP * UEqnP.H();
        surfaceScalarField phiHbyAP("phiHbyA", fvc::flux(HbyAP));
        adjustPhi(phiHbyAP, U, p);
        fvScalarMatrix pEqnP(fvm::laplacian(rAUP, p) == fvc::div(phiHbyAP));
        pEqnP.setReference(pRefCell, pRefValue);
        volScalarField pResP = pEqnP & p;
        surfaceScalarField phiResP = phiHbyAP - pEqnP.flux() - phi;

        // output the matrix-vector product to files
        word productName = "dRdWPsi_" + name(myProc) + "_FD_Values.txt";
        OFstream fOut(productName);
        forAll(URes, cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar val = (UResP[cellI][comp] - URes[cellI][comp]) / eps;
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
            scalar val = (pResP[cellI] - pRes[cellI]) / eps;
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
            scalar val = (phiResP[faceI] - phiRes[faceI]) / eps;
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
                scalar val = (phiResP.boundaryFieldRef()[patchI][faceI] - phiRes.boundaryFieldRef()[patchI][faceI]) / eps;
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
