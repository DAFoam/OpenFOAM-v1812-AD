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

Description
    Write primitive and binary block from OPstream

\*---------------------------------------------------------------------------*/

#include "UOPstream.H"
#include "PstreamGlobals.H"

#include <mpi.h>

// MediPack
#include <medi/medi.hpp>
#include <codi.hpp>
#include <codi/externals/codiMpiTypes.hpp>
using namespace medi;
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UOPstream::write
(
    const commsTypes commsType,
    const int toProcNo,
    const char* buf,
    const std::streamsize bufSize,
    const word callerInfo,
    const std::type_info& typeInfo,
    const int tag,
    const label communicator
)
{
    bool typeActive = Foam::UPstream::isTypeActive(typeInfo)
                   && codi::RealReverse::getGlobalTape().isActive();

    if (debug)
    {
        Pout<< "UOPstream::write : starting write to:" << toProcNo
            << " tag:" << tag
            << " comm:" << communicator << " size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;

        Pout<< " caller " << callerInfo
            << " typeActive: " << typeActive
            << " typeid: " << typeInfo.name()
            << Foam::endl;
    }
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "UOPstream::write : starting write to:" << toProcNo
            << " tag:" << tag
            << " comm:" << communicator << " size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Pout);
    }

    PstreamGlobals::checkCommunicator(communicator, toProcNo);

    bool transferFailed = true;
    // not checking the type
    if (commsType == commsTypes::blocking)
    {
        if(typeActive) 
        {
            transferFailed = AMPI_Bsend
            (
                reinterpret_cast<scalar*>(const_cast<char*>(buf)),
                bufSize/sizeof(scalar),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                toProcNo,   //procID(toProcNo),
                tag,
                PstreamGlobals::MPICommunicators_[communicator] //MPI_COMM_WORLD
            );
        } 
        else
        {
            transferFailed = AMPI_Bsend
            (
                const_cast<char*>(buf),
                bufSize,
                AMPI_BYTE,
                toProcNo,   //procID(toProcNo),
                tag,
                PstreamGlobals::MPICommunicators_[communicator] //MPI_COMM_WORLD
            );
        }

        if (debug)
        {
            Pout<< "UOPstream::write : finished write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }
    }
    else if (commsType == commsTypes::scheduled)
    {
        if(typeActive)
        {
            transferFailed = AMPI_Send
            (
                reinterpret_cast<scalar*>(const_cast<char*>(buf)),
                bufSize/sizeof(scalar),
                PstreamGlobals::mpiTypes_->MPI_TYPE, // MPI_BYTE,
                toProcNo,   //procID(toProcNo),
                tag,
                PstreamGlobals::MPICommunicators_[communicator] //MPI_COMM_WORLD
            );
        }
        else
        {
            transferFailed = AMPI_Send
            (
                const_cast<char*>(buf),
                bufSize,
                AMPI_BYTE, // AMPI_Type_No_Check, // MPI_BYTE,
                toProcNo,   //procID(toProcNo),
                tag,
                PstreamGlobals::MPICommunicators_[communicator] //MPI_COMM_WORLD
            );
        }

        if (debug)
        {
            Pout<< "UOPstream::write : finished write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }
    }
    else if (commsType == commsTypes::nonBlocking)
    {
        AMPI_Request request;
        if(typeActive)
        {
            transferFailed = AMPI_Isend
            (
                reinterpret_cast<scalar*>(const_cast<char*>(buf)),
                bufSize/sizeof(scalar),
                PstreamGlobals::mpiTypes_->MPI_TYPE, // MPI_BYTE,
                toProcNo,   //procID(toProcNo),
                tag,
                PstreamGlobals::MPICommunicators_[communicator],//MPI_COMM_WORLD,
                &request
            );
        }
        else
        {
            transferFailed = AMPI_Isend
            (
                const_cast<char*>(buf),
                bufSize,
                AMPI_BYTE,
                toProcNo,   //procID(toProcNo),
                tag,
                PstreamGlobals::MPICommunicators_[communicator],//MPI_COMM_WORLD,
                &request
            );
        }

        if (debug)
        {
            Pout<< "UOPstream::write : started write to:" << toProcNo
                << " tag:" << tag << " size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << " request:" << PstreamGlobals::outstandingRequests_.size()
                << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.append(request);
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << Foam::abort(FatalError);
    }

    return !transferFailed;
}

// ************************************************************************* //
