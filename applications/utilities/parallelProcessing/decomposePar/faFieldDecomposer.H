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
    faFieldDecomposer

Description
    Finite Area area and edge field decomposer.

Author
    Zeljko Tukovic, FSB Zagreb
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faFieldDecomposer.C
    faFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef faFieldDecomposer_H
#define faFieldDecomposer_H

#include "faMesh.H"
#include "faPatchFieldMapper.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class faFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class faFieldDecomposer
{
public:

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public faPatchFieldMapper
        {
            // Private data

                label sizeBeforeMapping_;
                labelList directAddressing_;

        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const label sizeBeforeMapping,
                    const labelUList& addressingSlice,
                    const label addressingOffset
                );


            // Member functions

                label size() const
                {
                    return directAddressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return true;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class
        class processorAreaPatchFieldDecomposer
        :
            public faPatchFieldMapper
        {
            // Private data

                label sizeBeforeMapping_;
                labelListList addressing_;
                scalarListList weights_;

        public:

            //- Construct given addressing
            processorAreaPatchFieldDecomposer
            (
                const faMesh& mesh,
                const labelUList& addressingSlice
            );


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return false;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


        //- Processor patch field decomposer class
        class processorEdgePatchFieldDecomposer
        :
            public faPatchFieldMapper
        {
            label sizeBeforeMapping_;
            labelListList addressing_;
            scalarListList weights_;

        public:

            //- Construct given addressing
            processorEdgePatchFieldDecomposer
            (
                label sizeBeforeMapping,
                const labelUList& addressingSlice
            );


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return false;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


private:

    // Private data

        //- Reference to complete mesh
        const faMesh& completeMesh_;

        //- Reference to processor mesh
        const faMesh& procMesh_;

        //- Reference to edge addressing
        const labelList& edgeAddressing_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        List<patchFieldDecomposer*> patchFieldDecomposerPtrs_;

        List<processorAreaPatchFieldDecomposer*>
            processorAreaPatchFieldDecomposerPtrs_;

        List<processorEdgePatchFieldDecomposer*>
            processorEdgePatchFieldDecomposerPtrs_;


    // Private Member Functions

        //- No copy construct
        faFieldDecomposer(const faFieldDecomposer&) = delete;

        //- No copy assignment
        void operator=(const faFieldDecomposer&) = delete;


public:

    // Constructors

        //- Construct from components
        faFieldDecomposer
        (
            const faMesh& completeMesh,
            const faMesh& procMesh,
            const labelList& edgeAddressing,
            const labelList& faceAddressing,
            const labelList& boundaryAddressing
        );


    // Destructor

        ~faFieldDecomposer();


    // Member Functions

        //- Decompose area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        decomposeField
        (
            const GeometricField<Type, faPatchField, areaMesh>& field
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        decomposeField
        (
            const GeometricField<Type, faePatchField, edgeMesh>& field
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
