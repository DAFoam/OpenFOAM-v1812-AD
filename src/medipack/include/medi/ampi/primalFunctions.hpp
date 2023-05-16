/*
 * MeDiPack, a Message Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of MeDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * MeDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * MeDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU
 * Lesser General Public License along with MeDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring (SciComp, University of Kaiserslautern-Landau)
 */

#pragma once

#include <iostream>

#include "ampiMisc.h"
#include "async.hpp"
#include "message.hpp"
#include "../displacementTools.hpp"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Send_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Send(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Isend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Isend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Bsend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Bsend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ibsend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Ibsend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ssend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Ssend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Issend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Issend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Rsend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Rsend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Irsend_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Irsend(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Recv_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int src, int tag, AMPI_Comm comm, AMPI_Status* status) {
    MEDI_UNUSED(count);
    MPI_Recv(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), src, tag, comm, status);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Mrecv_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, AMPI_Message* message, AMPI_Status* status) {
    MEDI_UNUSED(count);
    MPI_Mrecv(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), &message->message, status);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Irecv_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int src, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Irecv(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), src, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Imrecv_pri(typename DATATYPE::PrimalType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, AMPI_Message* message, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Imrecv(bufAdjoints, bufSize, datatype->getADTool().getPrimalMpiType(), &message->message, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Sendrecv_pri(typename SENDTYPE::PrimalType* sendbuf, int sendbufSize, int sendcount, SENDTYPE* sendtype, int dest, int sendtag,
                         typename RECVTYPE::PrimalType* recvbuf, int recvbufSize, int recvcount, RECVTYPE* recvtype, int source, int recvtag, AMPI_Comm comm, AMPI_Status*  status) {
#endif

    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);
    MPI_Sendrecv(sendbuf, sendbufSize, sendtype->getADTool().getPrimalMpiType(), dest, sendtag, recvbuf, recvbufSize, recvtype->getADTool().getPrimalMpiType(), source, recvtag, comm, status);
  }

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Bcast_wrap_pri(typename DATATYPE::PrimalType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(count);

    if(root == getCommRank(comm)) {
      std::swap(sendbufAdjoints, recvbufAdjoints);
    }
    MPI_Bcast(recvbufAdjoints, recvbufSize, datatype->getADTool().getPrimalMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_pri(typename DATATYPE::PrimalType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);

    if(root == getCommRank(comm)) {
      std::swap(sendbufAdjoints, recvbufAdjoints);
    }
    MPI_Ibcast(recvbufAdjoints, recvbufSize, datatype->getADTool().getPrimalMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatter_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Scatter(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Iscatter(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatterv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispl, MEDI_OPTIONAL_CONST int* sendcount, MEDI_OPTIONAL_CONST int* displs, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Scatterv(sendbufAdjoints, sendbufCounts, sendbufDispl, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispl, MEDI_OPTIONAL_CONST int* sendcount, MEDI_OPTIONAL_CONST int* displs, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Iscatterv(sendbufAdjoints, sendbufCounts, sendbufDispl, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gather_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Gather(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Igather(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gatherv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    MPI_Gatherv(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getPrimalMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    MPI_Igatherv(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getPrimalMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgather_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Allgather(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request ) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Iallgather(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgatherv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    MPI_Allgatherv(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getPrimalMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    MPI_Iallgatherv(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getPrimalMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoall_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Alltoall(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Ialltoall(sendbufAdjoints, sendbufSize, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, recvtype->getADTool().getPrimalMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoallv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispls, MEDI_OPTIONAL_CONST int* sendcounts, MEDI_OPTIONAL_CONST int* sdispls, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* rdispls, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcounts);
    MEDI_UNUSED(sdispls);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(rdispls);

    MPI_Alltoallv(sendbufAdjoints, sendbufCounts, sendbufDispls, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getPrimalMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_pri(typename SENDTYPE::PrimalType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispls, MEDI_OPTIONAL_CONST int* sendcounts, MEDI_OPTIONAL_CONST int* sdispls, SENDTYPE* sendtype, typename RECVTYPE::PrimalType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* rdispls, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcounts);
    MEDI_UNUSED(sdispls);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(rdispls);

    MPI_Ialltoallv(sendbufAdjoints, sendbufCounts, sendbufDispls, sendtype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getPrimalMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Reduce_global_pri(typename DATATYPE::PrimalType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, int root, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MEDI_UNUSED(sendbufSize);
    //TODO: MPI_Reduce(sendbufAdjoints, recvbufAdjoints, recvbufSize, datatype->getADTool().getPrimalMpiType(), TODO, root, comm);
    std::cout << "Forward reduce not supported." << std::endl;
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ireduce_global_pri(typename DATATYPE::PrimalType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MEDI_UNUSED(sendbufSize);
    //TODO: MPI_Ireduce(sendbufAdjoints, recvbufAdjoints, recvbufSize, datatype->getADTool().getPrimalMpiType(), TODO, root, comm, &request->request);
    std::cout << "Forward reduce not supported." << std::endl;
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Allreduce_global_pri(typename DATATYPE::PrimalType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MEDI_UNUSED(sendbufSize);

    //TODO MPI_Allreduce(sendbufAdjoints, datatype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, datatype->getADTool().getPrimalMpiType(), TODO, comm);
    std::cout << "Forward reduce not supported." << std::endl;
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Iallreduce_global_pri(typename DATATYPE::PrimalType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::PrimalType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(op);
    MEDI_UNUSED(count);

    //TODO MPI_Iallreduce(sendbufAdjoints, datatype->getADTool().getPrimalMpiType(), recvbufAdjoints, recvbufSize, datatype->getADTool().getPrimalMpiType(), comm, &request->request);
    std::cout << "Forward reduce not supported." << std::endl;
  }
#endif
}
