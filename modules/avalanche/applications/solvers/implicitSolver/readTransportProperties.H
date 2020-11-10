IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<frictionModel> friction
(
    frictionModel::New(transportProperties, Us, h, pb)
);

autoPtr<entrainmentModel> entrainment
(
    entrainmentModel::New(transportProperties, Us, h, hentrain, pb, tau)
);

autoPtr<depositionModel> deposition
(
    depositionModel::New(transportProperties, Us, h, hentrain, pb, tau)
);

dimensionedScalar hmin("hmin", dimLength, transportProperties);

dimensionedScalar xi("xi", dimless, transportProperties);

Switch pressureFeedback(transportProperties.lookup("pressureFeedback"));

Switch explicitDryAreas(transportProperties.lookup("explicitDryAreas"));

Switch bindHeight(transportProperties.lookupOrDefault<Switch>("bindHeight", true));

Info << "Running with" << endl
     << "    hmin " << hmin << endl
     << "    pressureFeedback is " << pressureFeedback << endl
     << "    explicitDryAreas is " << explicitDryAreas << endl
     << "    bindHeight is " << bindHeight << endl << endl;
