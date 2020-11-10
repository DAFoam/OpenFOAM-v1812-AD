/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::vtk::cloudAdaptor

Description
    A simple backend for converting an OpenFOAM cloud to vtkPolyData.

    The output is a multi-piece PolyData dataset corresponding to the
    cloud. Each piece corresponds to its MPI rank.

Note
    No caching of converted data since clouds are not stationary.

SourceFiles
    foamVtkCloudAdaptor.C
    foamVtkCloudAdaptorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkCloudAdaptor_H
#define foamVtkCloudAdaptor_H

#include "className.H"
#include "fvMesh.H"
#include "foamVtkTools.H"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkMultiPieceDataSet.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class vtk::cloudAdaptor Declaration
\*---------------------------------------------------------------------------*/

class cloudAdaptor
{
    // Private data

        const fvMesh& mesh_;


    // Private Member Functions

        //- Lagrangian positions as vtkPolyData
        static vtkSmartPointer<vtkPolyData> startLagrangian
        (
            const pointField& points
        );

        //- Lagrangian fields from objectRegistry
        template<class Type>
        static label convertLagrangianFields
        (
            vtkPolyData* vtkmesh,
            const objectRegistry& obr
        );

        //- Get cloud with point/cell data
        template<class UnaryPredicate>
        static vtkSmartPointer<vtkMultiPieceDataSet> getCloudImpl
        (
            const objectRegistry& mesh,
            const word& cloudName,
            const UnaryPredicate& pred
        );

    // Constructors

        //- No copy construct
        cloudAdaptor(const cloudAdaptor&) = delete;

        //- No copy assignment
        void operator=(const cloudAdaptor&) = delete;


public:

    //- Static Data Members

        ClassName("vtk::cloudAdaptor");


    // Constructors

        //- Construct for a particular mesh region
        cloudAdaptor(const fvMesh& mesh);


    //- Destructor
    ~cloudAdaptor() = default;


    // Member Functions

        //- Get cloud with point/cell data
        vtkSmartPointer<vtkMultiPieceDataSet> getCloud
        (
            const word& cloudName
        ) const;

        //- Get cloud with subset of point/cell data
        vtkSmartPointer<vtkMultiPieceDataSet> getCloud
        (
            const word& cloudName,
            const wordRes& selectFields
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
