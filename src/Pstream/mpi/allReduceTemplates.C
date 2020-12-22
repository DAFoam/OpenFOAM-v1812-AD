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

\*---------------------------------------------------------------------------*/

#include "allReduce.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, class BinaryOp>
void Foam::allReduce
(
    Type& Value,
    int MPICount,
    AMPI_Datatype MPIType,
    AMPI_Op AMPIOp,
    const BinaryOp& bop,
    const int tag,
    const label communicator
)
{
    if (!UPstream::parRun())
    {
        return;
    }
    
    if (UPstream::nProcs(communicator) <= UPstream::nProcsSimpleSum)
    {
        if (UPstream::master(communicator))
        {
            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave(communicator);
                slave++
            )
            {
                Type value;

                if
                (
                    AMPI_Recv
                    (
                        &value,
                        MPICount,
                        MPIType,
                        slave,  //UPstream::procID(slave),
                        tag,
                        PstreamGlobals::MPICommunicators_[communicator],
                        AMPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorInFunction
                        << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }

                Value = bop(Value, value);
            }
        }
        else
        {
            if
            (
                AMPI_Send
                (
                    &Value,
                    MPICount,
                    MPIType,
                    UPstream::masterNo(),//UPstream::procID(masterNo()),
                    tag,
                    PstreamGlobals::MPICommunicators_[communicator]
                )
            )
            {
                FatalErrorInFunction
                    << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }
        }


        if (UPstream::master(communicator))
        {
            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave(communicator);
                slave++
            )
            {
                if
                (
                    AMPI_Send
                    (
                        &Value,
                        MPICount,
                        MPIType,
                        slave,      //UPstream::procID(slave),
                        tag,
                        PstreamGlobals::MPICommunicators_[communicator]
                    )
                )
                {
                    FatalErrorInFunction
                        << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
            }
        }
        else
        {
            if
            (
                AMPI_Recv
                (
                    &Value,
                    MPICount,
                    MPIType,
                    UPstream::masterNo(),//UPstream::procID(masterNo()),
                    tag,
                    PstreamGlobals::MPICommunicators_[communicator],
                    AMPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorInFunction
                    << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
        }
    }
    else
    {
        Type sum;
        AMPI_Allreduce
        (
            &Value,
            &sum,
            MPICount,
            MPIType,
            AMPIOp,
            PstreamGlobals::MPICommunicators_[communicator]
        );
        Value = sum;
    }
}


// ************************************************************************* //
