/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "partialFaceAreaWeightAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
(
    label& startSeedi,
    label& srcFacei,
    label& tgtFacei,
    const boolList& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    const bool errorOnNotFound
) const
{
    faceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
    (
        startSeedi,
        srcFacei,
        tgtFacei,
        mapFlag,
        seedFaces,
        visitedFaces,
        false // no error on not found
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::
partialFaceAreaWeightAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
:
    faceAreaWeightAMI<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        triMode,
        reverseTarget,
        requireMatch
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::
~partialFaceAreaWeightAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
bool Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::conformal() const
{
    return false;
}


template<class SourcePatch, class TargetPatch>
void Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label srcFacei,
    label tgtFacei
)
{
    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            tgtAddress,
            tgtWeights,
            srcFacei,
            tgtFacei
        );

    if (!ok)
    {
        return;
    }

    // temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(this->srcPatch_.size());
    List<DynamicList<scalar>> srcWght(srcAddr.size());
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());
    List<DynamicList<scalar>> tgtWght(tgtAddr.size());

    faceAreaWeightAMI<SourcePatch, TargetPatch>::calcAddressing
    (
        srcAddr,
        srcWght,
        tgtAddr,
        tgtWght,
        srcFacei,
        tgtFacei
    );

    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i].transfer(srcWght[i]);
    }
    forAll(tgtAddr, i)
    {
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i].transfer(tgtWght[i]);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::setMagSf
(
    const TargetPatch& tgtPatch,
    const mapDistribute& map,
    scalarList& srcMagSf,
    scalarList& tgtMagSf
) const
{
    srcMagSf = std::move(this->srcMagSf_);

    scalarList newTgtMagSf(std::move(this->tgtMagSf_));
    map.reverseDistribute(tgtPatch.size(), newTgtMagSf);

    // Assign default sizes. Override selected values with
    // calculated values. This is to support ACMI
    // where some of the target faces are never used (so never get sent
    // over and hence never assigned to)
    tgtMagSf = tgtPatch.magFaceAreas();

    for (const labelList& smap : map.subMap())
    {
        UIndirectList<scalar>(tgtMagSf, smap) =
            UIndirectList<scalar>(newTgtMagSf, smap);
    }
}


// ************************************************************************* //
