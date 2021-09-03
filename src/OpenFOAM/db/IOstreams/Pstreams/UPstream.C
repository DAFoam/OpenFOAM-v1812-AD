/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "UPstream.H"
#include "debug.H"
#include "registerSwitch.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);
}

const Foam::Enum
<
    Foam::UPstream::commsTypes
>
Foam::UPstream::commsTypeNames
({
    { commsTypes::blocking, "blocking" },
    { commsTypes::scheduled, "scheduled" },
    { commsTypes::nonBlocking, "nonBlocking" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun(const label nProcs, const bool haveThreads)
{
    if (nProcs == 0)
    {
        parRun_ = false;
        haveThreads_ = haveThreads;

        freeCommunicator(UPstream::worldComm);
        label comm = allocateCommunicator(-1, labelList(1, label(0)), false);
        if (comm != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::worldComm:" << UPstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = "";
        Perr.prefix() = "";
    }
    else
    {
        parRun_ = true;
        haveThreads_ = haveThreads;

        // Redo worldComm communicator (this has been created at static
        // initialisation time)
        freeCommunicator(UPstream::worldComm);
        label comm = allocateCommunicator(-1, identity(nProcs), true);
        if (comm != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::worldComm:" << UPstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' +  name(myProcNo(Pstream::worldComm)) + "] ";
        Perr.prefix() = '[' +  name(myProcNo(Pstream::worldComm)) + "] ";
    }
}


Foam::label Foam::UPstream::allocateCommunicator
(
    const label parentIndex,
    const labelList& subRanks,
    const bool doPstream
)
{
    label index;
    if (!freeComms_.empty())
    {
        index = freeComms_.pop();
    }
    else
    {
        // Extend storage
        index = parentCommunicator_.size();

        myProcNo_.append(-1);
        procIDs_.append(List<int>(0));
        parentCommunicator_.append(-1);
        linearCommunication_.append(List<commsStruct>(0));
        treeCommunication_.append(List<commsStruct>(0));
    }

    if (debug)
    {
        Pout<< "Communicators : Allocating communicator " << index << endl
            << "    parent : " << parentIndex << endl
            << "    procs  : " << subRanks << endl
            << endl;
    }

    // Initialise; overwritten by allocatePstreamCommunicator
    myProcNo_[index] = 0;

    // Convert from label to int
    procIDs_[index].setSize(subRanks.size());
    forAll(procIDs_[index], i)
    {
        procIDs_[index][i] = subRanks[i];

        // Enforce incremental order (so index is rank in next communicator)
        if (i >= 1 && subRanks[i] <= subRanks[i-1])
        {
            FatalErrorInFunction
                << "subranks not sorted : " << subRanks
                << " when allocating subcommunicator from parent "
                << parentIndex
                << Foam::abort(FatalError);
        }
    }
    parentCommunicator_[index] = parentIndex;

    // Size but do not fill structure - this is done on-the-fly
    linearCommunication_[index] = List<commsStruct>(procIDs_[index].size());
    treeCommunication_[index] = List<commsStruct>(procIDs_[index].size());

    if (doPstream && parRun())
    {
        allocatePstreamCommunicator(parentIndex, index);
    }

    return index;
}


void Foam::UPstream::freeCommunicator
(
    const label communicator,
    const bool doPstream
)
{
    if (debug)
    {
        Pout<< "Communicators : Freeing communicator " << communicator << endl
            << "    parent   : " << parentCommunicator_[communicator] << endl
            << "    myProcNo : " << myProcNo_[communicator] << endl
            << endl;
    }

    if (doPstream && parRun())
    {
        freePstreamCommunicator(communicator);
    }
    myProcNo_[communicator] = -1;
    //procIDs_[communicator].clear();
    parentCommunicator_[communicator] = -1;
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    freeComms_.push(communicator);
}


void Foam::UPstream::freeCommunicators(const bool doPstream)
{
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freeCommunicator(communicator, doPstream);
        }
    }
}


int Foam::UPstream::baseProcNo(const label myComm, const int myProcID)
{
    int procID = myProcID;
    label comm = myComm;

    while (parent(comm) != -1)
    {
        const List<int>& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = UPstream::parent(comm);
    }

    return procID;
}


Foam::label Foam::UPstream::procNo(const label myComm, const int baseProcID)
{
    const List<int>& parentRanks = procID(myComm);
    label parentComm = parent(myComm);

    if (parentComm == -1)
    {
        return parentRanks.find(baseProcID);
    }
    else
    {
        const label parentRank = procNo(parentComm, baseProcID);
        return parentRanks.find(parentRank);
    }
}


Foam::label Foam::UPstream::procNo
(
    const label myComm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = UPstream::baseProcNo(currentComm, currentProcID);
    return procNo(myComm, physProcID);
}


template<>
Foam::UPstream::commsStruct&
Foam::UList<Foam::UPstream::commsStruct>::operator[](const label procID)
{
    UPstream::commsStruct& t = v_[procID];

    if (t.allBelow().size() + t.allNotBelow().size() + 1 != size())
    {
        // Not yet allocated

        label above(-1);
        labelList below(0);
        labelList allBelow(0);

        if (size() < UPstream::nProcsSimpleSum)
        {
            // Linear schedule

            if (procID == 0)
            {
                below.setSize(size()-1);
                for (label procI = 1; procI < size(); procI++)
                {
                    below[procI-1] = procI;
                }
            }
            else
            {
                above = 0;
            }
        }
        else
        {
            // Use tree like schedule. For 8 procs:
            // (level 0)
            //      0 receives from 1
            //      2 receives from 3
            //      4 receives from 5
            //      6 receives from 7
            // (level 1)
            //      0 receives from 2
            //      4 receives from 6
            // (level 2)
            //      0 receives from 4
            //
            // The sends/receives for all levels are collected per processor
            // (one send per processor; multiple receives possible) creating
            // a table:
            //
            // So per processor:
            // proc     receives from   sends to
            // ----     -------------   --------
            //  0       1,2,4           -
            //  1       -               0
            //  2       3               0
            //  3       -               2
            //  4       5               0
            //  5       -               4
            //  6       7               4
            //  7       -               6

            label mod = 0;

            for (label step = 1; step < size(); step = mod)
            {
                mod = step * 2;

                if (procID % mod)
                {
                    above = procID - (procID % mod);
                    break;
                }
                else
                {
                    for
                    (
                        label j = procID + step;
                        j < size() && j < procID + mod;
                        j += step
                    )
                    {
                        below.append(j);
                    }
                    for
                    (
                        label j = procID + step;
                        j < size() && j < procID + mod;
                        j++
                    )
                    {
                        allBelow.append(j);
                    }
                }
            }
        }
        t = UPstream::commsStruct(size(), procID, above, below, allBelow);
    }
    return t;
}


template<>
const Foam::UPstream::commsStruct&
Foam::UList<Foam::UPstream::commsStruct>::operator[](const label procID) const
{
    return const_cast<UList<UPstream::commsStruct>&>(*this).operator[](procID);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

Foam::LIFOStack<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<int> Foam::UPstream::myProcNo_(10);

Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(10);

Foam::DynamicList<Foam::label> Foam::UPstream::parentCommunicator_(10);

int Foam::UPstream::msgType_(1);


Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::linearCommunication_(10);

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::treeCommunication_(10);


// Allocate a serial communicator. This gets overwritten in parallel mode
// (by UPstream::setParRun())
Foam::UPstream::communicator serialComm
(
    -1,
    Foam::labelList(1, Foam::label(0)),
    false
);


bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 16)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.get
    (
        "commsType",
        Foam::debug::optimisationSwitches()
    )
);

// CoDiPack4OpenFOAM. The following variables are created to enable parallel AD
int Foam::UPstream::procOneToOneCommListIndex = -9999;

Foam::List< Foam::List<int> > Foam::UPstream::procOneToOneCommList;


void Foam::UPstream::calcOneToOneCommList
(
    List<DynamicList<label>> neighbProcList,
    List<List<label>>& commList
)
{
    /*
    Description:
        This function computes the processor one-to-one communication list.
        If we get a listList of neighbour processor indices, we want to compute
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
    // main iteration
    for(label iter=0;iter<maxIters;iter++)
    {
        if (iter==maxIters-1)
        {
            FatalErrorIn("calcOneToOneCommList")<<"maxIters reach! "<< exit(FatalError);
        }

        label finished = 1;
        List<label> empty;
        // append an empty list to commList for this iter
        commList.append(empty);
        //Info<<" commList "<<commList<<endl;
        // now loop over all processors to see if the procI is already in commList
        for(label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (commList[iter].found(procI))
            {
                // procI is already in commList, we need to check which neighbour procI is attached to this procI,
                // and we need to remove the neighbour procI from neighbProcList (it has been added in commList)
                // There are two possibilities:
                // 1. procI is an even number, which means it appears as neibProcI<->procI in commList
                // 2. procI is an odd number, which means it appears as procI<->neibProcI in commList
                // For case 1, we need to remove the index before procI, for case 2, we need to remove
                // the index after procI. This is because neibProcI<->procI always appear as a pair
                // See the following line in "else"
                //     commList[iter].append(procI);
                //     commList[iter].append(neighbProcI);
                label foundIndx = commList[iter].find(procI);
                label foundIndxRemainder = foundIndx % 2;
                label indxToRemove = -9999;
                if (foundIndxRemainder == 0)
                {
                    indxToRemove = foundIndx + 1;
                }
                else
                {
                    indxToRemove = foundIndx - 1;
                }
                label valToRemove = commList[iter][indxToRemove];
                forAll(neighbProcList[procI], idxI)
                {
                    const label& neiVal = neighbProcList[procI][idxI];
                    if (neiVal == valToRemove)
                    {
                        neighbProcList[procI].remove(idxI);
                        break;
                    }
                }
            }
            else
            {
                // procI is NOT in commList, we need to add procI into commList
                // we also need to add the neighbour procI that is associated with procI
                // into commList, if neighbour procI is not there yet.
                // Therefore, procI and neighbourProcI are always a pair in commList
                // After we add neighbour procI into commList, we need to remove 
                // it from  neighbProcList[procI]
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

namespace Foam
{
    // Register re-reader
    class addcommsTypeToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addcommsTypeToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addcommsTypeToOpt() = default;

        virtual void readData(Foam::Istream& is)
        {
            UPstream::defaultCommsType =
                UPstream::commsTypeNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os << UPstream::commsTypeNames[UPstream::defaultCommsType];
        }
    };

    addcommsTypeToOpt addcommsTypeToOpt_("commsType");
}

Foam::label Foam::UPstream::worldComm(0);

Foam::label Foam::UPstream::warnComm(-1);

int Foam::UPstream::nPollProcInterfaces
(
    Foam::debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitch
(
    "nPollProcInterfaces",
    int,
    Foam::UPstream::nPollProcInterfaces
);


int Foam::UPstream::maxCommsSize
(
    Foam::debug::optimisationSwitch("maxCommsSize", 0)
);
registerOptSwitch
(
    "maxCommsSize",
    int,
    Foam::UPstream::maxCommsSize
);


const int Foam::UPstream::mpiBufferSize
(
    Foam::debug::optimisationSwitch("mpiBufferSize", 0)
);


// ************************************************************************* //
