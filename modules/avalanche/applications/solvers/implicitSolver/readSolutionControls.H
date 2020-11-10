
    scalar nCorr =
            aMesh.solutionDict().lookupOrDefault<int>("nOuterCorrectors", 50);

    scalar minCorr =
            aMesh.solutionDict().lookupOrDefault<int>("minCorrectors", 3);

    scalar hResidualMax =
        aMesh.solutionDict().lookupOrDefault<scalar>("hResidualMax", 0.);

    scalar UsResidualMax =
        aMesh.solutionDict().lookupOrDefault<scalar>("UsResidualMax", 1);

