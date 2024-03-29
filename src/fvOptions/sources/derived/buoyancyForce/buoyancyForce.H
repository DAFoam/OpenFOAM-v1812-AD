/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::fv::buoyancyForce

Group
    grpFvOptionsSources

Description
    Calculates and applies the buoyancy force rho*g to the momentum equation
    corresponding to the specified velocity field.

Usage
    Example usage:
    \verbatim
    fields          (U);                    // Name of velocity field
    \endverbatim

SourceFiles
    buoyancyForce.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyForce_H
#define buoyancyForce_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyForce Declaration
\*---------------------------------------------------------------------------*/

class buoyancyForce
:
    public option
{
    // Private data

        uniformDimensionedVectorField g_;


    // Private Member Functions

        //- No copy construct
        buoyancyForce(const buoyancyForce&) = delete;

        //- No copy assignment
        void operator=(const buoyancyForce&) = delete;


public:

    //- Runtime type information
    TypeName("buoyancyForce");


    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Add explicit contribution to incompressible momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
