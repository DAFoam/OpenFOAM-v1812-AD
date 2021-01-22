/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Print caller information for debugging the MPI functions with AD
    Adapted from Towara 

\*---------------------------------------------------------------------------*/

#include <string>
#include <sstream>
#include <typeinfo>

inline std::string intToStr(int i){
    return std::to_string(i);
}

#ifdef __GNUC__
#define callerInfo() "in " + std::string(__FILE__) + ": " + std::string(__PRETTY_FUNCTION__) + std::string(":") + intToStr(__LINE__)
#else
#define callerInfo() "in " + std::string(__FILE__) + ": " + std::string(__FUNC__) + std::string(":") + intToStr(__LINE__)
#endif
