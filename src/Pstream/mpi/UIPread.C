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
    Read from UIPstream

\*---------------------------------------------------------------------------*/

#include "UIPstream.H"
#include "PstreamGlobals.H"
#include "IOstreams.H"

#include <mpi.h>

// MediPack
#include <medi/medi.hpp>
#include <codi.hpp>
#include <codi/externals/codiMpiTypes.hpp>
using namespace medi;
// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPstream::UIPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    DynamicList<char>& externalBuf,
    label& externalBufPosition,
    const int tag,
    const label comm,
    const bool clearAtEnd,
    streamFormat format,
    versionNumber version
)
:
    UPstream(commsType),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    externalBuf_(externalBuf),
    externalBufPosition_(externalBufPosition),
    tag_(tag),
    comm_(comm),
    clearAtEnd_(clearAtEnd),
    messageSize_(0)
{
    setOpened();
    setGood();

    if (commsType == commsTypes::nonBlocking)
    {
        // Message is already received into externalBuf
    }
    else
    {
        AMPI_Status status;

        label wantedSize = externalBuf_.capacity();

        if (debug)
        {
            Pout<< "UIPstream::UIPstream : read from:" << fromProcNo
                << " tag:" << tag << " comm:" << comm_
                << " wanted size:" << wantedSize
                << Foam::endl;
        }


        // If the buffer size is not specified, probe the incoming message
        // and set it
        if (!wantedSize)
        {
            AMPI_Probe
            (
                fromProcNo_,
                tag_,
                PstreamGlobals::MPICommunicators_[comm_],
                &status
            );
            AMPI_Get_count(&status, AMPI_BYTE, &messageSize_);

            externalBuf_.setCapacity(messageSize_);
            wantedSize = messageSize_;

            if (debug)
            {
                Pout<< "UIPstream::UIPstream : probed size:" << wantedSize
                    << Foam::endl;
            }
        }

        messageSize_ = UIPstream::read
        (
            commsType,
            fromProcNo_,
            externalBuf_.begin(),
            wantedSize,
            typeid(externalBuf_.begin()),
            tag_,
            comm_
        );

        // Set addressed size. Leave actual allocated memory intact.
        externalBuf_.setSize(messageSize_);

        if (!messageSize_)
        {
            setEof();
        }
    }
}


Foam::UIPstream::UIPstream(const int fromProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType_),
    Istream(buffers.format_, buffers.version_),
    fromProcNo_(fromProcNo),
    externalBuf_(buffers.recvBuf_[fromProcNo]),
    externalBufPosition_(buffers.recvBufPos_[fromProcNo]),
    tag_(buffers.tag_),
    comm_(buffers.comm_),
    clearAtEnd_(true),
    messageSize_(0)
{
    if
    (
        commsType() != UPstream::commsTypes::scheduled
     && !buffers.finishedSendsCalled_
    )
    {
        FatalErrorInFunction
            << "PstreamBuffers::finishedSends() never called." << endl
            << "Please call PstreamBuffers::finishedSends() after doing"
            << " all your sends (using UOPstream) and before doing any"
            << " receives (using UIPstream)" << Foam::exit(FatalError);
    }

    setOpened();
    setGood();

    if (commsType() == commsTypes::nonBlocking)
    {
        // Message is already received into externalBuf
        messageSize_ = buffers.recvBuf_[fromProcNo].size();

        if (debug)
        {
            Pout<< "UIPstream::UIPstream PstreamBuffers :"
                << " fromProcNo:" << fromProcNo
                << " tag:" << tag_ << " comm:" << comm_
                << " receive buffer size:" << messageSize_
                << Foam::endl;
        }
    }
    else
    {
        AMPI_Status status;

        label wantedSize = externalBuf_.capacity();

        if (debug)
        {
            Pout<< "UIPstream::UIPstream PstreamBuffers :"
                << " read from:" << fromProcNo
                << " tag:" << tag_ << " comm:" << comm_
                << " wanted size:" << wantedSize
                << Foam::endl;
        }

        // If the buffer size is not specified, probe the incoming message
        // and set it
        if (!wantedSize)
        {
            AMPI_Probe
            (
                fromProcNo_,
                tag_,
                PstreamGlobals::MPICommunicators_[comm_],
                &status
            );
            AMPI_Get_count(&status, AMPI_BYTE, &messageSize_);

            externalBuf_.setCapacity(messageSize_);
            wantedSize = messageSize_;

            if (debug)
            {
                Pout<< "UIPstream::UIPstream PstreamBuffers : probed size:"
                    << wantedSize << Foam::endl;
            }
        }

        messageSize_ = UIPstream::read
        (
            commsType(),
            fromProcNo_,
            externalBuf_.begin(),
            wantedSize,
            typeid(externalBuf_.begin()),
            tag_,
            comm_
        );

        // Set addressed size. Leave actual allocated memory intact.
        externalBuf_.setSize(messageSize_);

        if (!messageSize_)
        {
            setEof();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const std::type_info& typeInfo,
    const int tag,
    const label communicator
)
{

    bool typeActive = Foam::PstreamGlobals::isTypeActive(typeInfo)
                   && codi::RealReverse::getGlobalTape().isActive();

    if (debug)
    {
        Pout<< "UIPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << communicator
            << " wanted size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " typeActive: " << typeActive << " typeid: " << typeInfo.name()
            << Foam::endl;
    }
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "UIPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << communicator
            << " wanted size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Pout);
    }

    if (commsType == commsTypes::blocking || commsType == commsTypes::scheduled)
    {
        AMPI_Status status;
        label Err = 0;
        if(typeActive)
        {
            Err = AMPI_Recv
            (
                reinterpret_cast<scalar*>(buf),
                bufSize/sizeof(scalar),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            );
        }
        else
        {
            Err = AMPI_Recv
            (
                buf,
                bufSize,
                AMPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            );
        }
        if(Err)
        {
            FatalErrorInFunction
                << "MPI_Recv cannot receive incoming message"
                << Foam::abort(FatalError);

            return 0;
        }


        // Check size of message read

        int messageSize;
        AMPI_Get_count(&status, AMPI_BYTE, &messageSize);

        if (debug)
        {
            Pout<< "UIPstream::read : finished read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }

        if (messageSize > bufSize)
        {
            FatalErrorInFunction
                << "buffer (" << label(bufSize)
                << ") not large enough for incoming message ("
                << messageSize << ')'
                << Foam::abort(FatalError);
        }

        return messageSize;
    }
    else if (commsType == commsTypes::nonBlocking)
    {
        AMPI_Request request;
        label Err = 0;
        if (typeActive)
        {
            Err = AMPI_Irecv
            (
                reinterpret_cast<scalar*>(buf),
                bufSize/sizeof(scalar),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &request
            );
        }
        else
        {
            Err = AMPI_Irecv
            (
                buf,
                bufSize,
                AMPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &request
            );
        }

        if (Err)
        {
            FatalErrorInFunction
                << "MPI_Recv cannot start non-blocking receive"
                << Foam::abort(FatalError);

            return 0;
        }

        if (debug)
        {
            Pout<< "UIPstream::read : started read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << " request:" << PstreamGlobals::outstandingRequests_.size()
                << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.append(request);

        // Assume the message is completely received.
        return bufSize;
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << int(commsType)
            << Foam::abort(FatalError);

        return 0;
    }
}

// ************************************************************************* //
