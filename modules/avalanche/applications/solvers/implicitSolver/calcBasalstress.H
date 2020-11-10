
pb = -friction->rho()*xi*fac::ndiv(phi2s, Us)&n;

if (pressureFeedback)
{
    pb -= fac::ngrad((pbn&n)*h/2.)&n;
}

pb += friction->rho()*(g*h)&n;

//pb = fac::average(fac::interpolate(pb));

// Clipping
pb = Foam::max(pb, dimensionedScalar(dimPressure));

