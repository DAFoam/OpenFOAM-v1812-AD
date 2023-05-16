/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Test reverse mode partial derivatives dFdXv

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

    label myProc = Pstream::myProcNo();
    {
        codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
        tape.setActive();

        // compute dFdXv using reverse mode AD.

        // register input
        pointField meshPoints = mesh.points();
        forAll(meshPoints, pointI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                tape.registerInput(meshPoints[pointI][comp]);
            }
        }

        mesh.movePoints(meshPoints);

        volVectorField gradP = fvc::grad(p);
        scalar f = 0.0;

        if(Pstream::master())
        { 
            forAll(gradP, idxI) f+=mag(gradP[idxI]);
        }

        tape.registerOutput(f);
        // stop recording
        tape.setPassive();

        // Note: since we used reduced objFunc, we only need to
        // assign the seed for master proc
        //if (Pstream::master())
        //{
            f.setGradient(1.0);
        //}
        // evaluate tape to compute derivative
        tape.evaluate();

        // output the matrix-vector product to files
        word productName = "dFdXv_ProcI"  + name(myProc) + "_RAD.txt";
        OFstream fOut(productName);
        forAll(meshPoints, pointI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                scalar val = meshPoints[pointI][comp].getGradient();
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
