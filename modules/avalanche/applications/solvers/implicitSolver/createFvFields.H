volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity)
);

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength)
);

volScalarField Hentrain
(
    IOobject
    (
        "Hentrain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength)
);

volScalarField Pb
(
    IOobject
    (
        "Pb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimPressure)
);

volVectorField Tau
(
    IOobject
    (
        "Tau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimPressure)
);


#include "calcBasalstress.H"

// Create volume-to surface mapping object
volSurfaceMapping vsm(aMesh);
vsm.mapToVolume(Us, U.boundaryFieldRef());
vsm.mapToVolume(h, H.boundaryFieldRef());
vsm.mapToVolume(pb, Pb.boundaryFieldRef());
U.write();
H.write();
Pb.write();

