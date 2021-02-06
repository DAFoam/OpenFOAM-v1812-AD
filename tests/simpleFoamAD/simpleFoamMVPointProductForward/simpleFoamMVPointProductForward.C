/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Test forwad mode matrix vector product dRdXv * Psi and compare with FD

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
    argList::addOption(
        "mode",
        "FD",
        "Which mode to run");

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

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

    // AD
    if ("AD" == mode)
    {
        // compute dRdW * psi using forward mode AD. Here psi is a random vector
        // psi = sin(0.1*cellI) or sin(0.1*(cellI + comp))

        // set seeds
        word productNameSeed = "dRdXvPsi_" + name(myProc) + "_AD_Seeds.txt";
        OFstream fOutSeed(productNameSeed);

        pointField meshPoints = mesh.points();
        forAll(meshPoints, pointI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar randomSeed = sin(0.1 * (pointI + comp));
                fOutSeed << randomSeed << endl;
                meshPoints[pointI][comp].setGradient(randomSeed.getValue());
            }
        }
        mesh.movePoints(meshPoints);

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
        word productName = "dRdXvPsi_" + name(myProc) + "_AD_Values.txt";
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

    else if ("FD" == mode)
    {
        // compute dRdW * psi using finite-difference. Here psi is a random vector
        // psi = sin(0.1*cellI) or sin(0.1*(cellI + comp))

        // FD perturbation
        scalar eps = 1.0e-7;

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
        // perturb points
        pointField meshPoints = mesh.points();
        forAll(meshPoints, pointI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar randomSeed = sin(0.1 * (pointI + comp));
                meshPoints[pointI][comp] += randomSeed * eps;
            }
        }
        mesh.movePoints(meshPoints);

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
        word productName = "dRdXvPsi_" + name(myProc) + "_FD_Values.txt";
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

    Info<<"Done!"<<endl;
    return 0;
}

// ************************************************************************* //
