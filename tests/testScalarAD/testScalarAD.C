/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Compute wall distance for a given patch

\*---------------------------------------------------------------------------*/

#include "parRun.H"  
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar func(const scalar& a, const scalar& b)
{
    return 2.0*a+3.0*b*b;
}


int main(int argc, char *argv[])
{
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
     
    return 0;
}


// ************************************************************************* //
