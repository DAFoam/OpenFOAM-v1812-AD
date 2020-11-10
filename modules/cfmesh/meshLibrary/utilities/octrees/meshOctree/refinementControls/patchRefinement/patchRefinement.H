/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
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
    Foam::Module::patchRefinement

Description
    A class used for octree refinement. Refinement of boxes intersected
    by a given patch to the prescribed size.

SourceFiles
    patchRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef patchRefinement_H
#define patchRefinement_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class Istream;
class Ostream;

namespace Module
{
class triSurf;
class patchRefinement;
Istream& operator>>(Istream& is, patchRefinement& pr);
Ostream& operator<<(Ostream& os, const patchRefinement& pr);

/*---------------------------------------------------------------------------*\
                       Class patchRefinement Declaration
\*---------------------------------------------------------------------------*/

class patchRefinement
{
    // Private members

        word patchName_;

        scalar cellSize_;


public:

    // Constructors

        //- Null contructor
        patchRefinement();

        //- Construct from word and desired cell size
        patchRefinement(const word&, const scalar);

        //- contruct from stream
        patchRefinement(Istream&);


    //- Destructor
    ~patchRefinement() = default;


    // Member Functions

        const word& patchName() const;

        scalar cellSize() const;

        label patchInSurface(const triSurf& ts) const;


    // Operators

        void operator=(const patchRefinement&);

        friend Istream& operator>>(Istream& is, patchRefinement& pr);
        friend Ostream& operator<<(Ostream& os, const patchRefinement& pr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators
bool operator==(const patchRefinement& lhs, const patchRefinement& rhs);
bool operator!=(const patchRefinement& lhs, const patchRefinement& rhs);


} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
