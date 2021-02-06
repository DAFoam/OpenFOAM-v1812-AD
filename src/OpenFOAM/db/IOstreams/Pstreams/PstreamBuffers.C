/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "PstreamBuffers.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

Foam::DynamicList<char> Foam::PstreamBuffers::nullBuf(0);


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::PstreamBuffers
(
    const UPstream::commsTypes commsType,
    const word callerInfo,
    const bool typeActive,
    const int tag,
    const label comm,
    IOstream::streamFormat format,
    IOstream::versionNumber version
)
:
    commsType_(commsType),
    callerInfo_(callerInfo),
    typeActive_(typeActive),
    tag_(tag),
    comm_(comm),
    format_(format),
    version_(version),
    sendBuf_(UPstream::nProcs(comm)),
    recvBuf_(UPstream::nProcs(comm)),
    recvBufPos_(UPstream::nProcs(comm), 0),
    finishedSendsCalled_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    // Check that all data has been consumed.
    forAll(recvBufPos_, proci)
    {
        if (recvBufPos_[proci] < recvBuf_[proci].size())
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " not fully consumed. messageSize:" << recvBuf_[proci].size()
                << " bytes of which only " << recvBufPos_[proci]
                << " consumed."
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::finishedSends(const bool block)
{
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvBuf_,
            callerInfo_,
            typeActive_,
            tag_,
            comm_,
            block
        );
    }
}


void Foam::PstreamBuffers::finishedSends(labelList& recvSizes, const bool block)
{
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        Pstream::exchangeSizes(sendBuf_, recvSizes, comm_);

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvSizes,
            recvBuf_,
            callerInfo_,
            typeActive_,
            tag_,
            comm_,
            block
        );
    }
    else
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


void Foam::PstreamBuffers::clear()
{
    for (DynamicList<char>& buf : sendBuf_)
    {
        buf.clear();
    }
    for (DynamicList<char>& buf : recvBuf_)
    {
        buf.clear();
    }
    recvBufPos_ = 0;
    finishedSendsCalled_ = false;
}


void Foam::PstreamBuffers::calcOneToOneCommList
(
    List<DynamicList<label>> neighbProcList,
    List<List<label>>& commList
)
{
    /*
    Description:
        This function computes the one-to-one communication list. Essentially,
        if we get a listList of neighbour processor indices, we want to compute
        the commList such that, for each subList, the commList will have one
        processor communicating with only one other processor. This function is
        needed when exchanging information between processors. This is 
        because, for some reason, when we do nonBlocking comm, each processor can not 
        read/write data from/to more than one neighbour processor at a time.
        The message communication somehow mess up and we will get mixed (error) data.
        To fix this, we have to determine a so-called one-to-one communication
        list, such that, in this list, each subList contain the one-to-one
        processor communication index. Then we have to do commList.size() communication
        between processors, each time, we allow one processor to read/write data
        from only one neighbour and repeat this for all coupled patches.

    Example:
        If the neighbour processor list is like this:
        neighbProcList={
            {1,2,4},  // processor0's neighbour is 1, 2, and 4
            {0,3},    // processor1's neighbour is 0, and 3
            {0,3},
            {1,2,4},
            {0,3}
        };

        Then, the computed commList is (commList can't have repeated indices for each row)
        commList={
            {0,1,2,3}, // 1st round, we exchange data between procs 0<->1 and procs 2<->3
            {0,2,1,3}, // 2nd round, we exchange data between procs 0<->2 and procs 1<->3
            {0,4},
            {3,4}
        }
    */

    label maxIters = 1000;
    for(label iter=0;iter<maxIters;iter++)
    {
        if (iter==maxIters-1)
        {
            FatalErrorIn("calcOneToOneCommList")<<"maxIters reach! "<< exit(FatalError);
        }

        label finished = 1;
        List<label> empty;
        commList.append(empty);
        //Info<<" commList "<<commList<<endl;
        for(label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (commList[iter].found(procI))
            {
                forAll(neighbProcList[procI], idxI)
                {
                    const label& neighbProcI = neighbProcList[procI][idxI];
                    if(commList[iter].found(neighbProcI))
                    {
                        neighbProcList[procI].remove(idxI);
                        break;
                    }
                }
            }
            else
            {
 
                forAll(neighbProcList[procI], idxI)
                {
                    const label& neighbProcI = neighbProcList[procI][idxI];
                    if(!commList[iter].found(neighbProcI))
                    {
                        commList[iter].append(procI);
                        commList[iter].append(neighbProcI);
                        neighbProcList[procI].remove(idxI);
                        //Info<<" neighbProcList "<<neighbProcList<<endl;
                        //Info<<" commList "<<commList<<endl;
                        break;
                    }
                }

            }
            if (neighbProcList[procI].size()!=0)
            {
                finished=0;
            }
        }
        if(finished)
        {
            break;
        }
    }

}


// ************************************************************************* //
