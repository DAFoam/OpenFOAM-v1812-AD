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
        labelList testCellIs = {152, 182, 195};
        forAll(testCellIs, idxI)
        {
            label cellI = testCellIs[idxI];
            U[cellI][0].setGradient(1.0);

            fvVectorMatrix UEqn(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqn.relax();

            volVectorField URes = (UEqn & U) + fvc::grad(p);

            word partDerivResName = "verify_AD_CellI"+name(cellI);
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

        forAll(testCellIs, idxI)
        {
            label cellI = testCellIs[idxI];

            fvVectorMatrix UEqn(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqn.relax();

            volVectorField URes = (UEqn & U) + fvc::grad(p);

            U[cellI][0] += eps;
            U.correctBoundaryConditions();

            fvVectorMatrix UEqnP(
                fvm::div(phi, U)
                + turbulence->divDevReff(U));

            UEqnP.relax();

            volVectorField UResP = (UEqnP & U) + fvc::grad(p);

            U[cellI][0] -= eps;
            U.correctBoundaryConditions();
            
            word partDerivResName = "verify_FD_CellI"+name(cellI);
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
    else
    {
        // Test Proc0 CellI=61, 84
        // Proc1 CellI=59

    }

    return 0;
}


// ************************************************************************* //
