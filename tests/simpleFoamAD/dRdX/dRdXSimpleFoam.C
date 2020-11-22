/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Compute wall distance for a given patch

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "point",
        "0",
        "Which point to perturb"
    );
    argList::addOption
    (   
        "eps",
        "1e-4",
        "FD step"
    );
    argList::addOption
    (
        "proc",
        "0",
        "which procesor to perturb"
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    label myProc = Pstream::myProcNo();
    //label nProcs = Pstream::nProcs();

    label pointI = readLabel(args.optionLookup("point")());
    scalar eps = 1e-4;
    if (args.optionFound("eps"))
    {
        eps = readScalar(args.optionLookup("eps")());
    }
    label procI = 0;
    if (args.optionFound("proc"))
    {
        procI = readLabel(args.optionLookup("proc")());
    }
        // AD
        {
            pointField pointsNew = mesh.points();
            if(procI == myProc)
            {
                pointsNew[pointI][1].setGradient(1.0);
            }
            mesh.movePoints(pointsNew);

            // Res
            fvVectorMatrix UEqn(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqn.relax();

            volVectorField URes = (UEqn & U) + fvc::grad(p);

            word partDerivResName = "verify_ProcI_"+name(myProc)+"_PointI_"+name(pointI)+"_AD.txt";
            OFstream fOut(partDerivResName);
            label glbIdx = 0;
            forAll(URes, idxJ)
            {
                for(label j=0;j<3;j++)
                {
                    fOut<<"glbIdx: "<<glbIdx;
                    glbIdx++;
                    scalar deriv = URes[idxJ][j].getGradient();
                    if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                    else fOut<<" deriv: 0";
                    fOut<<" var: U"<<j;
                    fOut<<" x: "<<mesh.C()[idxJ][0]<<" y: "<<mesh.C()[idxJ][1]<<" z: "<<mesh.C()[idxJ][2]<<endl;
                }
            }

            label pRefCell = 0;
            scalar pRefValue = 0.0;

            volScalarField rAU(1.0 / UEqn.A());
            volVectorField HbyA("HbyA", U);
            HbyA = rAU * UEqn.H();

            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

            adjustPhi(phiHbyA, U, p);

            fvScalarMatrix pEqn(
                fvm::laplacian(rAU, p)
                == fvc::div(phiHbyA));
            pEqn.setReference(pRefCell, pRefValue);

            volScalarField pRes = pEqn & p;
            surfaceScalarField phiRes = phiHbyA - pEqn.flux() - phi;

            forAll(pRes, idxJ)
            {
                fOut<<"glbIdx: "<<glbIdx;
                glbIdx++;
                scalar deriv = pRes[idxJ].getGradient();
                if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                else fOut<<" deriv: 0";
                fOut<<" var: p";
                fOut<<" x: "<<mesh.C()[idxJ][0]<<" y: "<<mesh.C()[idxJ][1]<<" z: "<<mesh.C()[idxJ][2]<<endl;
            }
            forAll(phiRes, idxJ)
            {       
                fOut<<"glbIdx: "<<glbIdx;
                glbIdx++;
                scalar deriv = phiRes[idxJ].getGradient();
                if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                else fOut<<" deriv: 0";
                fOut<<" var: phi";
                fOut<<" x: "<<mesh.Cf()[idxJ][0]<<" y: "<<mesh.Cf()[idxJ][1]<<" z: "<<mesh.Cf()[idxJ][2]<<endl;
            }

        }

        // FD
        {
            // ref Res
            fvVectorMatrix UEqn(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqn.relax();

            volVectorField URes = (UEqn & U) + fvc::grad(p);

            label pRefCell = 0;
            scalar pRefValue = 0.0;

            volScalarField rAU(1.0 / UEqn.A());
            volVectorField HbyA("HbyA", U);
            HbyA = rAU * UEqn.H();

            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

            adjustPhi(phiHbyA, U, p);

            fvScalarMatrix pEqn(
                fvm::laplacian(rAU, p)
                == fvc::div(phiHbyA));
            pEqn.setReference(pRefCell, pRefValue);

            volScalarField pRes = pEqn & p;
            surfaceScalarField phiRes = phiHbyA - pEqn.flux() - phi;

            // perturb
            pointField pointsNew = mesh.points();
            if (procI == myProc) {
                pointsNew[pointI][1] += eps;
            }
            mesh.movePoints(pointsNew);

            // perturbed res
            fvVectorMatrix UEqnP(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqnP.relax();

            volVectorField UResP = (UEqnP & U) + fvc::grad(p);

            volScalarField rAUP(1.0 / UEqnP.A());
            volVectorField HbyAP("HbyA", U);
            HbyAP = rAUP * UEqnP.H();

            surfaceScalarField phiHbyAP("phiHbyA", fvc::flux(HbyAP));

            adjustPhi(phiHbyAP, U, p);

            fvScalarMatrix pEqnP(
                fvm::laplacian(rAUP, p)
                == fvc::div(phiHbyAP));
            pEqnP.setReference(pRefCell, pRefValue);

            volScalarField pResP = pEqnP & p;
            surfaceScalarField phiResP = phiHbyAP - pEqnP.flux() - phi;
            
            word partDerivResName = "verify_ProcI_"+name(myProc)+"_PointI_"+name(pointI)+"_FD.txt";
            OFstream fOut(partDerivResName);
            label glbIdx = 0;
            forAll(URes, idxJ)
            {
                for(label j=0;j<3;j++)
                {
                    fOut<<"glbIdx: "<<glbIdx;
                    glbIdx++;
                    scalar deriv = (UResP[idxJ][j] - URes[idxJ][j] ) /eps;
                    if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                    else fOut<<" deriv: 0";
                    fOut<<" var: U"<<j;
                    fOut<<" x: "<<mesh.C()[idxJ][0]<<" y: "<<mesh.C()[idxJ][1]<<" z: "<<mesh.C()[idxJ][2]<<endl;
                }
            }

            forAll(pRes, idxJ)
            {   
                fOut<<"glbIdx: "<<glbIdx;
                glbIdx++;
                scalar deriv = (pResP[idxJ] - pRes[idxJ] ) /eps;
                if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                else fOut<<" deriv: 0";
                fOut<<" var: p";
                fOut<<" x: "<<mesh.C()[idxJ][0]<<" y: "<<mesh.C()[idxJ][1]<<" z: "<<mesh.C()[idxJ][2]<<endl;
            }
            forAll(phiRes, idxJ)
            {       
                fOut<<"glbIdx: "<<glbIdx;
                glbIdx++;
                scalar deriv = (phiResP[idxJ] - phiRes[idxJ] ) /eps;
                if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                else fOut<<" deriv: 0";
                fOut<<" var: phi";
                fOut<<" x: "<<mesh.Cf()[idxJ][0]<<" y: "<<mesh.Cf()[idxJ][1]<<" z: "<<mesh.Cf()[idxJ][2]<<endl;
            }    

        }     


    return 0;
}


// ************************************************************************* //
