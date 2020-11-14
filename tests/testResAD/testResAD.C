/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Compute wall distance for a given patch

\*---------------------------------------------------------------------------*/

#include "fvCFD.H" 
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	/*
    scalar a = 1.0;
    scalar b = 2.0;

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();
    tape.registerInput(a);
    tape.registerInput(b);
    scalar c = func(a,b);
    tape.registerOutput(c);
    tape.setPassive();
    c.setGradient(1.0);
    tape.evaluate();

    Info<<c.getValue()<<endl;
    Info<<a.getGradient()<<endl;
    Info<<b.getGradient()<<endl;
    */
     #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    return 0;
}


// ************************************************************************* //
