/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::doubleScalar

Description
    Floating-point double precision scalar type.

SourceFiles
    doubleScalar.C

\*---------------------------------------------------------------------------*/

#ifndef doubleScalar_H
#define doubleScalar_H

#include "doubleFloat.H"
#include "direction.H"
#include "word.H"
// CoDiPack4OpenFOAM  Add CoDiPack header
#include "codi.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// CoDiPack4OpenFOAM  Replace double with codi::RealReverse
typedef codi::RealReverse doubleScalar;

// Largest and smallest scalar values allowed in certain parts of the code.
// (15 is the number of significant figures in an
//  IEEE double precision number.  See limits.h or float.h)
// CoDiPack4OpenFOAM  Change doubleScalar to double
constexpr double doubleScalarGREAT = 1.0e+15;
constexpr double doubleScalarVGREAT = 1.0e+300;
constexpr double doubleScalarROOTVGREAT = 1.0e+150;
constexpr double doubleScalarSMALL = 1.0e-15;
constexpr double doubleScalarROOTSMALL = 3.0e-8;
constexpr double doubleScalarVSMALL = 1.0e-300;
constexpr double doubleScalarROOTVSMALL = 1.0e-150;


#define Scalar doubleScalar
#define ScalarVGREAT doubleScalarVGREAT
#define ScalarVSMALL doubleScalarVSMALL
#define ScalarROOTVGREAT doubleScalarROOTVGREAT
#define ScalarROOTVSMALL doubleScalarROOTVSMALL
#define ScalarRead readDouble


inline Scalar mag(const Scalar s)
{
    // CoDiPack4OpenFOAM
    return fabs(s);
}

inline Scalar hypot(const Scalar x, const Scalar y)
{
    // CoDiPack4OpenFOAM  Manually define hypot function 
    // Ref: http://en.cppreference.com/w/cpp/numeric/math/hypot
    return sqrt(x*x+y*y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    // CoDiPack4OpenFOAM
    return atan2(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    // CoDiPack4OpenFOAM
    return jn(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    // CoDiPack4OpenFOAM
    return yn(n, s);
}

// Normal (double-precision) transcendental functions
// CoDiPack4OpenFOAM
#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return func(s);                \
}


// Double-precision bessel functions
// CoDiPack4OpenFOAM
#define besselFunc(func)           \
inline Scalar func(const Scalar s) \
{                                  \
    return func(s);                \
}

// CoDiPack4OpenFOAM Explicit instantiation of operator<< for double
Ostream& operator<<(Ostream& os, const double s);
Istream& operator>>(Istream& is, double& val);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scalar.H"

#undef Scalar
#undef ScalarVGREAT
#undef ScalarVSMALL
#undef ScalarROOTVGREAT
#undef ScalarROOTVSMALL
#undef ScalarRead
#undef transFunc
#undef besselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
