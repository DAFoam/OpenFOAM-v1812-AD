/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

Class
    Foam::volSurfaceMapping

Description
    Volume to surface and surface to volume mapping

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    volSurfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef volSurfaceMapping_H
#define volSurfaceMapping_H

#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvPatchField;

/*---------------------------------------------------------------------------*\
                      Class volSurfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class volSurfaceMapping
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        volSurfaceMapping(const volSurfaceMapping&) = delete;

        //- No copy assignment
        void operator=(const volSurfaceMapping&) = delete;


public:

    // Constructors

        //- Construct null
        volSurfaceMapping();

        //- Construct from mesh
        volSurfaceMapping(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Destructor - default


    // Member Functions

        //- Map droplet cloud sources to surface
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const typename
            GeometricField<Type, fvPatchField, volMesh>::Boundary& df
        ) const;

        //- Map surface field to volume boundary field
        template<class Type>
        void mapToVolume
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bf
        ) const;

        template<class Type>
        void mapToVolume
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volSurfaceMapping.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
