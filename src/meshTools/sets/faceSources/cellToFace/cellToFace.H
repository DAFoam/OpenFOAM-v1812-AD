/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cellToFace

Description
    A topoSetFaceSource to select a faceSet from a cellSet.

    Either use 'all' cell faces, or only faces that have cells
    on 'both' sides.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        sets        | The cellSet names to use          | possibly |
        set         | The cellSet name to use           | possibly |
        option      | Selection (all/both)              | yes   |
    \endtable

Note
    Must specify "sets" or "set" (highest to lowest precedence).
    When picking up coupled faces, it uses cells on neighbouring processors.

SourceFiles
    cellToFace.C

\*---------------------------------------------------------------------------*/

#ifndef cellToFace_H
#define cellToFace_H

#include "topoSetFaceSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellToFace Declaration
\*---------------------------------------------------------------------------*/

class cellToFace
:
    public topoSetFaceSource
{
public:
        //- Enumeration defining the valid options
        enum cellAction
        {
            ALL,
            BOTH
        };


private:

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<cellAction> cellActionNames_;

        //- Names of sets to use
        wordList names_;

        //- Option
        cellAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("cellToFace");

    // Constructors

        //- Construct from components
        cellToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const cellAction option
        );

        //- Construct from dictionary
        cellToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cellToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cellToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
