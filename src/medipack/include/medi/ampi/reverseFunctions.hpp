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
  void AMPI_Send_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Recv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, MPI_STATUS_IGNORE);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Isend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Irecv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Bsend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Recv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, MPI_STATUS_IGNORE);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ibsend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Irecv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ssend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Recv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, MPI_STATUS_IGNORE);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Issend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Irecv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Rsend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    MPI_Recv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, MPI_STATUS_IGNORE);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Irsend_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int dest, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Irecv(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), dest, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Recv_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int src, int tag, AMPI_Comm comm, AMPI_Status* status) {
    MEDI_UNUSED(count);
    MEDI_UNUSED(status);
    MPI_Send(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), src, tag, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Mrecv_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, AMPI_Message* message, AMPI_Status* status) {
    MEDI_UNUSED(count);
    MEDI_UNUSED(status);
    MPI_Send(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), message->src, message->tag, message->comm);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Irecv_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, int src, int tag, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Isend(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), src, tag, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Imrecv_adj(typename DATATYPE::AdjointType* bufAdjoints, int bufSize, int count, DATATYPE* datatype, AMPI_Message* message, AMPI_Request* request) {
    MEDI_UNUSED(count);
    MPI_Isend(bufAdjoints, bufSize, datatype->getADTool().getAdjointMpiType(), message->src, message->tag, message->comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Sendrecv_adj(typename SENDTYPE::AdjointType* sendbuf, int sendbufSize, int sendcount, SENDTYPE* sendtype, int dest, int sendtag,
#endif
                     typename RECVTYPE::AdjointType* recvbuf, int recvbufSize, int recvcount, RECVTYPE* recvtype, int source, int recvtag, AMPI_Comm comm, AMPI_Status*  status) {

    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);
    MPI_Sendrecv(recvbuf, recvbufSize, recvtype->getADTool().getAdjointMpiType(), source, recvtag, sendbuf, sendbufSize, sendtype->getADTool().getAdjointMpiType(), dest, sendtag, comm, status);
  }

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Bcast_wrap_adj(typename DATATYPE::AdjointType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(count);

    MPI_Gather(recvbufAdjoints, recvbufSize, datatype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, datatype->getADTool().getAdjointMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_adj(typename DATATYPE::AdjointType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);

    MPI_Igather(recvbufAdjoints, recvbufSize, datatype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, datatype->getADTool().getAdjointMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatter_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Gather(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Igather(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatterv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispl, MEDI_OPTIONAL_CONST int* sendcount, MEDI_OPTIONAL_CONST int* displs, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Gatherv(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufCounts, sendbufDispl, sendtype->getADTool().getAdjointMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispl, MEDI_OPTIONAL_CONST int* sendcount, MEDI_OPTIONAL_CONST int* displs, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Igatherv(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufCounts, sendbufDispl, sendtype->getADTool().getAdjointMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gather_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Scatter(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Iscatter(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gatherv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    MPI_Scatterv(recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), root, comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    MPI_Iscatterv(recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), root, comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgather_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Alltoall(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request ) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Ialltoall(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgatherv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    LinearDisplacements linDis(getCommSize(comm), sendbufSize);

    MPI_Alltoallv(recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, linDis.counts, linDis.displs, sendtype->getADTool().getAdjointMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(displs);

    LinearDisplacements* linDis = new LinearDisplacements(getCommSize(comm), sendbufSize);
    request->setReverseData(reinterpret_cast<void*>(linDis), LinearDisplacements::deleteFunc);

    MPI_Ialltoallv(recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, linDis->counts, linDis->displs, sendtype->getADTool().getAdjointMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoall_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Alltoall(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int sendbufSize, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcount);
    MEDI_UNUSED(recvcount);

    MPI_Ialltoall(recvbufAdjoints, recvbufSize, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, sendtype->getADTool().getAdjointMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoallv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispls, MEDI_OPTIONAL_CONST int* sendcounts, MEDI_OPTIONAL_CONST int* sdispls, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* rdispls, RECVTYPE* recvtype, AMPI_Comm comm) {
    MEDI_UNUSED(sendcounts);
    MEDI_UNUSED(sdispls);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(rdispls);

    MPI_Alltoallv(recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufCounts, sendbufDispls, sendtype->getADTool().getAdjointMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_adj(typename SENDTYPE::AdjointType* &sendbufAdjoints, int* sendbufCounts, MEDI_OPTIONAL_CONST int* sendbufDispls, MEDI_OPTIONAL_CONST int* sendcounts, MEDI_OPTIONAL_CONST int* sdispls, SENDTYPE* sendtype, typename RECVTYPE::AdjointType* &recvbufAdjoints, int* recvbufCounts, MEDI_OPTIONAL_CONST int* recvbufDispls, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* rdispls, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(sendcounts);
    MEDI_UNUSED(sdispls);
    MEDI_UNUSED(recvcounts);
    MEDI_UNUSED(rdispls);

    MPI_Ialltoallv(recvbufAdjoints, recvbufCounts, recvbufDispls, recvtype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufCounts, sendbufDispls, sendtype->getADTool().getAdjointMpiType(), comm, &request->request);
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Reduce_global_adj(typename DATATYPE::AdjointType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, int root, AMPI_Comm comm) {
    MEDI_UNUSED(count);
    if(root == getCommRank(comm)) {
      MPI_Bcast(recvbufAdjoints, recvbufSize, datatype->getADTool().getAdjointMpiType(), root, comm);
      std::swap(sendbufAdjoints, recvbufAdjoints);
    } else {
      MPI_Bcast(sendbufAdjoints, sendbufSize, datatype->getADTool().getAdjointMpiType(), root, comm);
    }
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Ireduce_global_adj(typename DATATYPE::AdjointType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, int root, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(count);
    if(root == getCommRank(comm)) {
      MPI_Ibcast(recvbufAdjoints, recvbufSize, datatype->getADTool().getAdjointMpiType(), root, comm, &request->request);
      std::swap(sendbufAdjoints, recvbufAdjoints);
    } else {
      MPI_Ibcast(sendbufAdjoints, sendbufSize, datatype->getADTool().getAdjointMpiType(), root, comm, &request->request);
    }
  }
#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Allreduce_global_adj(typename DATATYPE::AdjointType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, AMPI_Comm comm) {
    MEDI_UNUSED(op);
    MEDI_UNUSED(count);

    MPI_Allgather(recvbufAdjoints, recvbufSize, datatype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, datatype->getADTool().getAdjointMpiType(), comm);
  }
#endif

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  void AMPI_Iallreduce_global_adj(typename DATATYPE::AdjointType* &sendbufAdjoints, int sendbufSize, typename DATATYPE::AdjointType* &recvbufAdjoints, int recvbufSize, int count, DATATYPE* datatype, AMPI_Op op, AMPI_Comm comm, AMPI_Request* request) {
    MEDI_UNUSED(op);
    MEDI_UNUSED(count);

    MPI_Iallgather(recvbufAdjoints, recvbufSize, datatype->getADTool().getAdjointMpiType(), sendbufAdjoints, sendbufSize, datatype->getADTool().getAdjointMpiType(), comm, &request->request);
  }
#endif
}
