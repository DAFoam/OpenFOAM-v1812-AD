/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Typedef
    Foam::scalar

Description
    A floating-point number identical to float or double depending on
    whether WM_SP or WM_DP is defined.

SourceFiles
    scalar.C

\*---------------------------------------------------------------------------*/

#ifndef scalar_H
#define scalar_H

#include "floatScalar.H"
#include "doubleScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP)

// Define scalar as a float

namespace Foam
{
    typedef floatScalar scalar;
    static const scalar GREAT = floatScalarGREAT;
    static const scalar VGREAT = floatScalarVGREAT;
    static const scalar ROOTVGREAT = floatScalarROOTVGREAT;
    static const scalar SMALL = floatScalarSMALL;
    static const scalar ROOTSMALL = floatScalarROOTSMALL;
    static const scalar VSMALL = floatScalarVSMALL;
    static const scalar ROOTVSMALL = floatScalarROOTVSMALL;

    scalar readScalar(Istream& is);

    inline scalar readScalar(const char* buf)
    {
        return readFloat(buf);
    }

    inline bool readScalar(const char* buf, scalar& val)
    {
        return readFloat(buf, val);
    }

    inline scalar readScalar(const std::string& str)
    {
        return readFloat(str);
    }

    inline bool readScalar(const std::string& str, scalar& val)
    {
        return readFloat(str, val);
    }
}

#elif defined(WM_DP)

// Define scalar as a double

namespace Foam
{
    typedef doubleScalar scalar;
    static const scalar GREAT = doubleScalarGREAT;
    static const scalar VGREAT = doubleScalarVGREAT;
    static const scalar ROOTVGREAT = doubleScalarROOTVGREAT;
    static const scalar SMALL = doubleScalarSMALL;
    static const scalar ROOTSMALL = doubleScalarROOTSMALL;
    static const scalar VSMALL = doubleScalarVSMALL;
    static const scalar ROOTVSMALL = doubleScalarROOTVSMALL;

    scalar readScalar(Istream& is);

    inline scalar readScalar(const char* buf)
    {
        return readDouble(buf);
    }

    inline bool readScalar(const char* buf, scalar& val)
    {
        return readDouble(buf, val);
    }

    inline scalar readScalar(const std::string& str)
    {
        return readDouble(str);
    }

    inline bool readScalar(const std::string& str, scalar& val)
    {
        return readDouble(str, val);
    }
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Additional transcendental functions

namespace Foam
{
    //- Inverse normalized incomplete gamma function
    scalar invIncGamma(const scalar a, const scalar P);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
