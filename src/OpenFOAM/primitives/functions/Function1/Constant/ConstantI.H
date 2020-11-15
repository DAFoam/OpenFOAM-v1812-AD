/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "Constant.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1Types::Constant<Type>::value(const scalar x) const
{
    return value_;
}

// CoDiPack4OpenFOAM check this
// all return types except label
template <typename R, typename T>
typename std::enable_if<!std::is_same<R, Foam::label>::value, R>::type
    passive_if_label_scalar(const T& t) {
        return t;
}
// if return type label
template <typename R, typename T>
typename std::enable_if<std::is_same<R, Foam::label>::value, R>::type
    passive_if_label_scalar(const T& t) {
        return t.getValue();
}


template<class Type>
inline Type Foam::Function1Types::Constant<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // CoDiPack4OpenFOAM check this
    return passive_if_label_scalar<Type>((x2 - x1)*value_);
}


// ************************************************************************* //
