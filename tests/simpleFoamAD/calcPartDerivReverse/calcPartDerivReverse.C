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
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    label myProc = Pstream::myProcNo();
    label nProcs = Pstream::nProcs();

    scalar eps = 1e-6;

    if (nProcs == 1)
    {
        // Test CellI 152 172 195
        labelList testCellIs = {152}; //, 182, 195};
        forAll(testCellIs, idxI)
        {

            codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
            tape.setActive();

            forAll(U, idxJ)
            {   
                for(label j=0;j<3;j++) tape.registerInput(U[idxJ][j]);
                tape.registerInput(p[idxJ]);
            }
            forAll(phi, idxJ) tape.registerInput(phi[idxJ]);

            // Res
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


            label cellI = testCellIs[idxI];
            tape.registerOutput(URes[cellI][0]);
            tape.setPassive();
            URes[cellI][0].setGradient(1.0);
            tape.evaluate();

            word partDerivResName = "verify_AD_CellI"+name(cellI);
            OFstream fOut(partDerivResName);
            label glbIdx = 0;
            forAll(U, idxJ)
            {
                for(label j=0;j<3;j++)
                {
                    fOut<<"glbIdx: "<<glbIdx;
                    glbIdx++;
                    scalar deriv = U[idxJ][j].getGradient();
                    if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                    else fOut<<" deriv: 0";
                    fOut<<" var: U"<<j;
                    fOut<<" x: "<<mesh.C()[idxJ][0]<<" y: "<<mesh.C()[idxJ][1]<<" z: "<<mesh.C()[idxJ][2]<<endl;
                }
            }

            forAll(p, idxJ)
            {
                fOut<<"glbIdx: "<<glbIdx;
                glbIdx++;
                scalar deriv = p[idxJ].getGradient();
                if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                else fOut<<" deriv: 0";
                fOut<<" var: p";
                fOut<<" x: "<<mesh.C()[idxJ][0]<<" y: "<<mesh.C()[idxJ][1]<<" z: "<<mesh.C()[idxJ][2]<<endl;
            }
            forAll(phi, idxJ)
            {       
                fOut<<"glbIdx: "<<glbIdx;
                glbIdx++;
                scalar deriv = phi[idxJ].getGradient();
                if (fabs(deriv) > 1e-16) fOut<<" deriv: "<<deriv;
                else fOut<<" deriv: 0";
                fOut<<" var: phi";
                fOut<<" x: "<<mesh.Cf()[idxJ][0]<<" y: "<<mesh.Cf()[idxJ][1]<<" z: "<<mesh.Cf()[idxJ][2]<<endl;
            }

        }

        forAll(testCellIs, idxI)
        {
            label cellI = testCellIs[idxI];

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


            U[cellI][0] += eps;
            U.correctBoundaryConditions();

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

            U[cellI][0] -= eps;
            U.correctBoundaryConditions();
            
            word partDerivResName = "verify_FD_CellI"+name(cellI);
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

    }
    else
    {
        // Test Proc0 CellI=61, 84
        // Proc1 CellI=59

        // cellI = 61
        label cellI = 61;
        {
            if (myProc == 0)
            {
                U[cellI][0].setGradient(1.0);
            }
    
            fvVectorMatrix UEqn(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));
            UEqn.relax();
    
            volVectorField URes = (UEqn & U) + fvc::grad(p);
    
            word partDerivResName = "verify_AD_CellI"+name(cellI)+"_n"+name(myProc);
            OFstream fOut(partDerivResName);
            forAll(URes, idxJ)
            {
                for(label j=0;j<3;j++)
                {
                    scalar deriv = URes[idxJ][j].getGradient();
                    if (fabs(deriv) > 1e-16) fOut<<deriv<<endl;
                    else fOut<<"0"<<endl;
                }
            }
        }

        {
            fvVectorMatrix UEqn(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));
            UEqn.relax();
    
            volVectorField URes = (UEqn & U) + fvc::grad(p);

            if (myProc == 0)
            {
                U[cellI][0]+=eps;
            }
            U.correctBoundaryConditions();

            fvVectorMatrix UEqnP(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqnP.relax();

            volVectorField UResP = (UEqnP & U) + fvc::grad(p);

            if (myProc == 0)
            {
                U[cellI][0] -= eps;
            }
            U.correctBoundaryConditions();
            
            word partDerivResName = "verify_FD_CellI"+name(cellI)+"_n"+name(myProc);
            OFstream fOut(partDerivResName);
            forAll(URes, idxJ)
            {
                for(label j=0;j<3;j++)
                {
                    scalar deriv = (UResP[idxJ][j] - URes[idxJ][j] ) /eps;
                    if (fabs(deriv) > 1e-16) fOut<<deriv<<endl;
                    else fOut<<"0"<<endl;
                }
            }    
        }

    }

    return 0;
}


// ************************************************************************* //
