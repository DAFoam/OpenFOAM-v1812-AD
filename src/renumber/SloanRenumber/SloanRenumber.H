/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

Class
    Foam::SloanRenumber

Description
    Sloan renumbering algorithm

    E.g. http://www.apav.it/sito_ratio/file_pdf/ratio_4/capitolo_1.pdf

SourceFiles
    SloanRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef SloanRenumber_H
#define SloanRenumber_H

#include "renumberMethod.H"
#include "Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SloanRenumber Declaration
\*---------------------------------------------------------------------------*/

class SloanRenumber
:
    public renumberMethod
{
    // Private data

        const bool reverse_;


    // Private Member Functions

        //- No copy construct
        SloanRenumber(const SloanRenumber&) = delete;

        //- No copy assignment
        void operator=(const SloanRenumber&) = delete;


public:

    //- Runtime type information
    TypeName("Sloan");


    // Constructors

        //- Construct given the renumber dictionary
        SloanRenumber(const dictionary& renumberDict);


    //- Destructor
    virtual ~SloanRenumber() = default;


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
            return labelList(0);
        }

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
