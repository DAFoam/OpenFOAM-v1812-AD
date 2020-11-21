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
     #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

	label resI = 100;
	scalar resRef = 0.0;
	
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

	forAll(p, cellI)
{
                p[cellI] = cos(1.0*cellI);
}
p.correctBoundaryConditions();

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

forAll(U, cellI)
{
	for(label i=0;i<3;i++)
	{
		U[cellI][i] = sin(1.0*cellI+i);
}
}
U.correctBoundaryConditions();

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);

        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();

    tape.registerInput(U[100][0]);

fvVectorMatrix UEqn
    (
     fvm::div(phi, U) - fvc::grad(p)
    );

volVectorField URes = UEqn & U;
   
resRef = URes[resI][0];
tape.registerOutput(URes[resI][0]);
    tape.setPassive();
    URes[resI][0].setGradient(1.0);
    tape.evaluate();

                   Info<<"AD: "<< U[100][0].getGradient()<<endl;

scalar eps =1e-8;
    U[100][0] += eps;

    fvVectorMatrix UEqnP
    (
     fvm::div(phi, U) - fvc::grad(p)
    );

volVectorField UResP = UEqnP & U;

scalar resP = UResP[resI][0];

Info<<"FD: "<< (resP-resRef)/eps<<endl;

    return 0;
}


// ************************************************************************* //
