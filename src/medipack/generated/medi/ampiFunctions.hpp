/*
 * MeDiPack, a Message Differentiation Package
 *
 * Copyright (C) 2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum (SciComp, TU Kaiserslautern)
 *
 * This file is part of MeDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * MeDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * MeDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with MeDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "../../include/medi/ampi/async.hpp"
#include "../../include/medi/ampi/message.hpp"
#include "../../include/medi/ampi/reverseFunctions.hpp"
#include "../../include/medi/ampi/forwardFunctions.hpp"
#include "../../include/medi/ampi/primalFunctions.hpp"
#include "../../include/medi/displacementTools.hpp"
#include "../../include/medi/mpiTools.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Bsend_AdjointHandle : public HandleBase {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Bsend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Bsend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Bsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Bsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Bsend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Bsend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Bsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Bsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Bsend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Bsend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Bsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Bsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Bsend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Bsend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                 AMPI_Comm comm) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Bsend(buf, count, datatype->getMpiType(), dest, tag, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Bsend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Bsend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Bsend_b<DATATYPE>;
        h->funcForward = AMPI_Bsend_d<DATATYPE>;
        h->funcPrimal = AMPI_Bsend_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Bsend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm);
      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Ibsend_AdjointHandle : public AsyncAdjointHandle {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Ibsend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Ibsend_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Ibsend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Ibsend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ibsend_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Ibsend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Ibsend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ibsend_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Ibsend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Ibsend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ibsend_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Ibsend_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Ibsend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                  AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ibsend(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ibsend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Ibsend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Ibsend_b<DATATYPE>;
        h->funcForward = AMPI_Ibsend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Ibsend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Ibsend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Ibsend_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Ibsend_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Ibsend_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Ibsend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Ibsend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Ibsend_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ibsend_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Ibsend_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET

  template<typename DATATYPE>
  struct AMPI_Bsend_init_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };


  template<typename DATATYPE>
  int AMPI_Bsend_init_preStart(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Bsend_init_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Bsend_init_postEnd(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Bsend_init(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                      AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Bsend_init(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ibsend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      rStatus = MPI_Bsend_init(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Bsend_init_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Bsend_init_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Bsend_init_finish<DATATYPE>;
      request->start = (ContinueFunction)AMPI_Bsend_init_preStart<DATATYPE>;
      request->end = (ContinueFunction)AMPI_Bsend_init_postEnd<DATATYPE>;
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Bsend_init_preStart(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Bsend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Bsend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      int bufElements = 0;

      // recompute the total size of the buffer
      bufElements = count;
      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Ibsend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Ibsend_b<DATATYPE>;
        h->funcForward = AMPI_Ibsend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Ibsend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      asyncHandle->toolHandle = h;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Ibsend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Ibsend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Bsend_init_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Bsend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Bsend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Bsend_init_postEnd(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Bsend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Bsend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ibsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {


      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Imrecv_AdjointHandle : public AsyncAdjointHandle {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    typename DATATYPE::PrimalType* bufOldPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    AMPI_Message message;

    ~AMPI_Imrecv_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
      if(nullptr != bufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufOldPrimals);
        bufOldPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Imrecv_AsyncHandle : public AsyncHandle {
    typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    AMPI_Message* message;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Imrecv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );

    AMPI_Imrecv_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, &h->message, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Imrecv_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Imrecv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Imrecv_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, &h->message, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Imrecv_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Imrecv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }

    AMPI_Imrecv_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, &h->message, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Imrecv_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Imrecv_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Imrecv(typename DATATYPE::Type* buf, int count, DATATYPE* datatype, AMPI_Message* message,
                  AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Imrecv(buf, count, datatype->getMpiType(), &message->message, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Imrecv_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Imrecv_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufOldPrimals, h->bufTotalSize);
          datatype->getValues(buf, 0, h->bufOldPrimals, 0, count);
        }



        datatype->createIndices(buf, 0, h->bufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Imrecv_b<DATATYPE>;
        h->funcForward = AMPI_Imrecv_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Imrecv_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->message = *message;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(buf, 0, count);
      }

      rStatus = MPI_Imrecv(bufMod, count, datatype->getModifiedMpiType(), &message->message, &request->request);

      AMPI_Imrecv_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Imrecv_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->message = message;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Imrecv_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Imrecv_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Imrecv_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Imrecv_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Imrecv_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Imrecv_AsyncHandle<DATATYPE>*>(handle);
    typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    AMPI_Message* message = asyncHandle->message;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Imrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Imrecv_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(message); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(buf, 0, h->bufIndices, h->bufOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Irecv_AdjointHandle : public AsyncAdjointHandle {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    typename DATATYPE::PrimalType* bufOldPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int source;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Irecv_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
      if(nullptr != bufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufOldPrimals);
        bufOldPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Irecv_AsyncHandle : public AsyncHandle {
    typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int source;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Irecv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );

    AMPI_Irecv_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->source, h->tag, h->comm,
                             &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Irecv_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Irecv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Irecv_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->source, h->tag, h->comm,
                             &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Irecv_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Irecv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }

    AMPI_Irecv_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->source, h->tag, h->comm,
                             &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Irecv_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Irecv_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Irecv(typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int source, int tag, AMPI_Comm comm,
                 AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Irecv(buf, count, datatype->getMpiType(), source, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Irecv_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Irecv_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufOldPrimals, h->bufTotalSize);
          datatype->getValues(buf, 0, h->bufOldPrimals, 0, count);
        }



        datatype->createIndices(buf, 0, h->bufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Irecv_b<DATATYPE>;
        h->funcForward = AMPI_Irecv_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Irecv_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->source = source;
        h->tag = tag;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(buf, 0, count);
      }

      rStatus = MPI_Irecv(bufMod, count, datatype->getModifiedMpiType(), source, tag, comm, &request->request);

      AMPI_Irecv_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Irecv_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->source = source;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Irecv_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Irecv_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Irecv_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Irecv_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Irecv_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Irecv_AsyncHandle<DATATYPE>*>(handle);
    typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int source = asyncHandle->source;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(source); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(buf, 0, h->bufIndices, h->bufOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Irsend_AdjointHandle : public AsyncAdjointHandle {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Irsend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Irsend_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Irsend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Irsend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Irsend_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Irsend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Irsend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Irsend_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Irsend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Irsend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Irsend_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Irsend_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Irsend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                  AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Irsend(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Irsend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Irsend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Irsend_b<DATATYPE>;
        h->funcForward = AMPI_Irsend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Irsend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Irsend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Irsend_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Irsend_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Irsend_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Irsend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Irsend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Irsend_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Irsend_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Irsend_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Isend_AdjointHandle : public AsyncAdjointHandle {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Isend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Isend_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Isend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Isend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                             &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Isend_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Isend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Isend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                             &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Isend_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Isend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Isend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                             &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Isend_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Isend_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Isend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                 AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Isend(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Isend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Isend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Isend_b<DATATYPE>;
        h->funcForward = AMPI_Isend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Isend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Isend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Isend_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Isend_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Isend_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Isend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Isend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Isend_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Isend_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Isend_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Issend_AdjointHandle : public AsyncAdjointHandle {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Issend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Issend_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Issend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Issend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Issend_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Issend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Issend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Issend_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Issend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Issend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm,
                              &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Issend_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Issend_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Issend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                  AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Issend(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Issend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Issend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Issend_b<DATATYPE>;
        h->funcForward = AMPI_Issend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Issend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Issend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Issend_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Issend_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Issend_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Issend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Issend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Issend_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Issend_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Issend_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Mrecv_AdjointHandle : public HandleBase {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    typename DATATYPE::PrimalType* bufOldPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    AMPI_Message message;
    AMPI_Status* status;

    ~AMPI_Mrecv_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
      if(nullptr != bufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufOldPrimals);
        bufOldPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Mrecv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Mrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Mrecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );

    AMPI_Mrecv_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, &h->message, h->status);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Mrecv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Mrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Mrecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Mrecv_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, &h->message, h->status);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Mrecv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Mrecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Mrecv_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }

    AMPI_Mrecv_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, &h->message, h->status);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Mrecv(typename DATATYPE::Type* buf, int count, DATATYPE* datatype, AMPI_Message* message,
                 AMPI_Status* status) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Mrecv(buf, count, datatype->getMpiType(), &message->message, status);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Mrecv_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Mrecv_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufOldPrimals, h->bufTotalSize);
          datatype->getValues(buf, 0, h->bufOldPrimals, 0, count);
        }



        datatype->createIndices(buf, 0, h->bufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Mrecv_b<DATATYPE>;
        h->funcForward = AMPI_Mrecv_d<DATATYPE>;
        h->funcPrimal = AMPI_Mrecv_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->message = *message;
        h->status = status;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(buf, 0, count);
      }

      rStatus = MPI_Mrecv(bufMod, count, datatype->getModifiedMpiType(), &message->message, status);
      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(buf, 0, h->bufIndices, h->bufOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Recv_AdjointHandle : public HandleBase {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    typename DATATYPE::PrimalType* bufOldPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int source;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Recv_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
      if(nullptr != bufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufOldPrimals);
        bufOldPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Recv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Recv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Recv_AdjointHandle<DATATYPE>*>(handle);

    MPI_Status status;
    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );

    AMPI_Recv_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->source, h->tag, h->comm, &status);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Recv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Recv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Recv_AdjointHandle<DATATYPE>*>(handle);

    MPI_Status status;
    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Recv_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->source, h->tag, h->comm, &status);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Recv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Recv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Recv_AdjointHandle<DATATYPE>*>(handle);

    MPI_Status status;
    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->bufIndices, h->bufOldPrimals, h->bufTotalSize);
    }

    AMPI_Recv_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->source, h->tag, h->comm, &status);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Recv(typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int source, int tag, AMPI_Comm comm,
                AMPI_Status* status) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Recv(buf, count, datatype->getMpiType(), source, tag, comm, status);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Recv_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Recv_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufOldPrimals, h->bufTotalSize);
          datatype->getValues(buf, 0, h->bufOldPrimals, 0, count);
        }



        datatype->createIndices(buf, 0, h->bufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Recv_b<DATATYPE>;
        h->funcForward = AMPI_Recv_d<DATATYPE>;
        h->funcPrimal = AMPI_Recv_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->source = source;
        h->tag = tag;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(buf, 0, count);
      }

      rStatus = MPI_Recv(bufMod, count, datatype->getModifiedMpiType(), source, tag, comm, status);
      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(buf, 0, h->bufIndices, h->bufOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET

  template<typename DATATYPE>
  struct AMPI_Recv_init_AsyncHandle : public AsyncHandle {
    typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int source;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };


  template<typename DATATYPE>
  int AMPI_Recv_init_preStart(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Recv_init_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Recv_init_postEnd(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Recv_init(typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int source, int tag, AMPI_Comm comm,
                     AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Recv_init(buf, count, datatype->getMpiType(), source, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Irecv_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      rStatus = MPI_Recv_init(bufMod, count, datatype->getModifiedMpiType(), source, tag, comm, &request->request);

      AMPI_Recv_init_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Recv_init_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->source = source;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Recv_init_finish<DATATYPE>;
      request->start = (ContinueFunction)AMPI_Recv_init_preStart<DATATYPE>;
      request->end = (ContinueFunction)AMPI_Recv_init_postEnd<DATATYPE>;
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Recv_init_preStart(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Recv_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Recv_init_AsyncHandle<DATATYPE>*>(handle);
    typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int source = asyncHandle->source;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(source); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      int bufElements = 0;

      // recompute the total size of the buffer
      bufElements = count;
      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Irecv_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufOldPrimals, h->bufTotalSize);
          datatype->getValues(buf, 0, h->bufOldPrimals, 0, count);
        }



        datatype->createIndices(buf, 0, h->bufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Irecv_b<DATATYPE>;
        h->funcForward = AMPI_Irecv_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Irecv_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->source = source;
        h->tag = tag;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(buf, 0, count);
      }

      asyncHandle->toolHandle = h;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Irecv_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Irecv_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Recv_init_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Recv_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Recv_init_AsyncHandle<DATATYPE>*>(handle);
    typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int source = asyncHandle->source;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(source); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(buf, 0, h->bufIndices, h->bufOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Recv_init_postEnd(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Recv_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Recv_init_AsyncHandle<DATATYPE>*>(handle);
    typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int source = asyncHandle->source;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irecv_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irecv_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(source); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {


      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Rsend_AdjointHandle : public HandleBase {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Rsend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Rsend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Rsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Rsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Rsend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Rsend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Rsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Rsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Rsend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Rsend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Rsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Rsend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Rsend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Rsend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                 AMPI_Comm comm) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Rsend(buf, count, datatype->getMpiType(), dest, tag, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Rsend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Rsend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Rsend_b<DATATYPE>;
        h->funcForward = AMPI_Rsend_d<DATATYPE>;
        h->funcPrimal = AMPI_Rsend_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Rsend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm);
      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET

  template<typename DATATYPE>
  struct AMPI_Rsend_init_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };


  template<typename DATATYPE>
  int AMPI_Rsend_init_preStart(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Rsend_init_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Rsend_init_postEnd(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Rsend_init(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                      AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Rsend_init(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Irsend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      rStatus = MPI_Rsend_init(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Rsend_init_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Rsend_init_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Rsend_init_finish<DATATYPE>;
      request->start = (ContinueFunction)AMPI_Rsend_init_preStart<DATATYPE>;
      request->end = (ContinueFunction)AMPI_Rsend_init_postEnd<DATATYPE>;
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Rsend_init_preStart(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Rsend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Rsend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      int bufElements = 0;

      // recompute the total size of the buffer
      bufElements = count;
      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Irsend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Irsend_b<DATATYPE>;
        h->funcForward = AMPI_Irsend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Irsend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      asyncHandle->toolHandle = h;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Irsend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Irsend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Rsend_init_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Rsend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Rsend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Rsend_init_postEnd(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Rsend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Rsend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Irsend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Irsend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {


      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Send_AdjointHandle : public HandleBase {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Send_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Send_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Send_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Send_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Send_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Send_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Send_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Send_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Send_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Send_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Send_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Send_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Send_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Send(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                AMPI_Comm comm) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Send(buf, count, datatype->getMpiType(), dest, tag, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Send_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Send_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Send_b<DATATYPE>;
        h->funcForward = AMPI_Send_d<DATATYPE>;
        h->funcPrimal = AMPI_Send_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Send(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm);
      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET

  template<typename DATATYPE>
  struct AMPI_Send_init_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };


  template<typename DATATYPE>
  int AMPI_Send_init_preStart(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Send_init_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Send_init_postEnd(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Send_init(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                     AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Send_init(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Isend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      rStatus = MPI_Send_init(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Send_init_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Send_init_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Send_init_finish<DATATYPE>;
      request->start = (ContinueFunction)AMPI_Send_init_preStart<DATATYPE>;
      request->end = (ContinueFunction)AMPI_Send_init_postEnd<DATATYPE>;
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Send_init_preStart(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Send_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Send_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      int bufElements = 0;

      // recompute the total size of the buffer
      bufElements = count;
      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Isend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Isend_b<DATATYPE>;
        h->funcForward = AMPI_Isend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Isend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      asyncHandle->toolHandle = h;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Isend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Isend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Send_init_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Send_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Send_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Send_init_postEnd(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Send_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Send_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Isend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Isend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {


      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Sendrecv_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int dest;
    int sendtag;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int source;
    int recvtag;
    AMPI_Comm comm;

    ~AMPI_Sendrecv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Sendrecv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    MPI_Status status;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Sendrecv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype, h->dest,
                                          h->sendtag, h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->source, h->recvtag, h->comm, &status);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Sendrecv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    MPI_Status status;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Sendrecv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype, h->dest,
                                          h->sendtag, h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->source, h->recvtag, h->comm, &status);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Sendrecv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    MPI_Status status;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Sendrecv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype, h->dest,
                                          h->sendtag, h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->source, h->recvtag, h->comm, &status);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Sendrecv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype, int dest,
                    int sendtag, typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, int source, int recvtag,
                    AMPI_Comm comm, AMPI_Status* status) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Sendrecv(sendbuf, sendcount, sendtype->getMpiType(), dest, sendtag, recvbuf, recvcount,
                             recvtype->getMpiType(), source, recvtag, comm, status);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      sendbufElements = sendcount;

      if(sendtype->isModifiedBufferRequired() ) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = recvcount;

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Sendrecv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->sendbufCount = sendtype->computeActiveElements(sendcount);
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = recvtype->computeActiveElements(recvcount);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount);
        }


        sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);

        recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Sendrecv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Sendrecv_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Sendrecv_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->dest = dest;
        h->sendtag = sendtag;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->source = source;
        h->recvtag = recvtag;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        recvtype->clearIndices(recvbuf, 0, recvcount);
      }

      rStatus = MPI_Sendrecv(sendbufMod, sendcount, sendtype->getModifiedMpiType(), dest, sendtag, recvbufMod, recvcount,
                             recvtype->getModifiedMpiType(), source, recvtag, comm, status);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount);
      }

      if(nullptr != h) {
        // handle the recv buffers
        recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount);
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired() ) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Ssend_AdjointHandle : public HandleBase {
    int bufTotalSize;
    typename DATATYPE::IndexType* bufIndices;
    typename DATATYPE::PrimalType* bufPrimals;
    /* required for async */ void* bufAdjoints;
    int bufCount;
    int bufCountVec;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;

    ~AMPI_Ssend_AdjointHandle () {
      if(nullptr != bufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufIndices);
        bufIndices = nullptr;
      }
      if(nullptr != bufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufPrimals);
        bufPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Ssend_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ssend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ssend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufPrimals, h->bufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->bufIndices, h->bufPrimals, h->bufTotalSize);


    AMPI_Ssend_pri<DATATYPE>(h->bufPrimals, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Ssend_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ssend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ssend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);


    AMPI_Ssend_fwd<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Ssend_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ssend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ssend_AdjointHandle<DATATYPE>*>(handle);

    h->bufAdjoints = nullptr;
    h->bufCountVec = adjointInterface->getVectorSize() * h->bufCount;
    adjointInterface->createAdjointTypeBuffer(h->bufAdjoints, h->bufTotalSize );

    AMPI_Ssend_adj<DATATYPE>(h->bufAdjoints, h->bufCountVec, h->count, h->datatype, h->dest, h->tag, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufIndices, h->bufAdjoints, h->bufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Ssend(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                 AMPI_Comm comm) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ssend(buf, count, datatype->getMpiType(), dest, tag, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ssend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Ssend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Ssend_b<DATATYPE>;
        h->funcForward = AMPI_Ssend_d<DATATYPE>;
        h->funcPrimal = AMPI_Ssend_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      rStatus = MPI_Ssend(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm);
      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET

  template<typename DATATYPE>
  struct AMPI_Ssend_init_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf;
    typename DATATYPE::ModifiedType* bufMod;
    int count;
    DATATYPE* datatype;
    int dest;
    int tag;
    AMPI_Comm comm;
    AMPI_Request* request;
  };


  template<typename DATATYPE>
  int AMPI_Ssend_init_preStart(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Ssend_init_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Ssend_init_postEnd(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Ssend_init(MEDI_OPTIONAL_CONST typename DATATYPE::Type* buf, int count, DATATYPE* datatype, int dest, int tag,
                      AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ssend_init(buf, count, datatype->getMpiType(), dest, tag, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Issend_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufMod = nullptr;
      int bufElements = 0;

      // compute the total size of the buffer
      bufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufMod, bufElements);
      } else {
        bufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(buf));
      }

      rStatus = MPI_Ssend_init(bufMod, count, datatype->getModifiedMpiType(), dest, tag, comm, &request->request);

      AMPI_Ssend_init_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Ssend_init_AsyncHandle<DATATYPE>();
      asyncHandle->buf = buf;
      asyncHandle->bufMod = bufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->dest = dest;
      asyncHandle->tag = tag;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Ssend_init_finish<DATATYPE>;
      request->start = (ContinueFunction)AMPI_Ssend_init_preStart<DATATYPE>;
      request->end = (ContinueFunction)AMPI_Ssend_init_postEnd<DATATYPE>;
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Ssend_init_preStart(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ssend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Ssend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      int bufElements = 0;

      // recompute the total size of the buffer
      bufElements = count;
      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Issend_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        datatype->copyIntoModifiedBuffer(buf, 0, bufMod, 0, count);
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        h->bufCount = datatype->computeActiveElements(count);
        h->bufTotalSize = datatype->computeActiveElements(bufElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufIndices, h->bufTotalSize);




        datatype->getIndices(buf, 0, h->bufIndices, 0, count);


        // pack all the variables in the handle
        h->funcReverse = AMPI_Issend_b<DATATYPE>;
        h->funcForward = AMPI_Issend_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Issend_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->dest = dest;
        h->tag = tag;
        h->comm = comm;
      }


      asyncHandle->toolHandle = h;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Issend_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Issend_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Ssend_init_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ssend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Ssend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings


    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);


      if(nullptr != h) {
        // handle the recv buffers
      }

      datatype->getADTool().stopAssembly(h);
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Ssend_init_postEnd(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ssend_init_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Ssend_init_AsyncHandle<DATATYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* buf = asyncHandle->buf;
    typename DATATYPE::ModifiedType* bufMod = asyncHandle->bufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int dest = asyncHandle->dest;
    int tag = asyncHandle->tag;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Issend_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Issend_AdjointHandle<DATATYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(buf); // Unused generated to ignore warnings
    MEDI_UNUSED(bufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(dest); // Unused generated to ignore warnings
    MEDI_UNUSED(tag); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {


      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Allgather_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Allgather_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgather_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Allgather_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                           h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgather_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Allgather_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                           h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgather_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize * getCommSize(h->comm));

    AMPI_Allgather_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                           h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm);

    adjointInterface->combineAdjoints(h->sendbufAdjoints, h->sendbufTotalSize, getCommSize(h->comm));
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Allgather(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                     typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Allgather(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcount;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = recvcount * getCommSize(comm);

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Allgather_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          recvtype->copyIntoModifiedBuffer(recvbuf, recvcount * getCommRank(comm), recvbufMod, recvcount * getCommRank(comm),
                                           recvcount);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(recvcount);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = recvtype->computeActiveElements(recvcount);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          recvtype->getIndices(recvbuf, recvcount * getCommRank(comm), h->sendbufIndices, 0, recvcount);
        }

        recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount * getCommSize(comm));

        // pack all the variables in the handle
        h->funcReverse = AMPI_Allgather_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Allgather_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Allgather_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        recvtype->clearIndices(recvbuf, 0, recvcount * getCommSize(comm));
      }

      rStatus = MPI_Allgather(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                              recvtype->getModifiedMpiType(), comm);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
      }

      if(nullptr != h) {
        // handle the recv buffers
        recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Allgatherv_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int* recvbufCount;
    /* required for async */ int* recvbufCountVec;
    /* required for async */ int* recvbufDisplsVec;
    MEDI_OPTIONAL_CONST  int* recvcounts;
    MEDI_OPTIONAL_CONST  int* displs;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Allgatherv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
      if(nullptr != recvbufCount) {
        delete [] recvbufCount;
        recvbufCount = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgatherv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Allgatherv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                            h->recvbufPrimals, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgatherv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Allgatherv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                            h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Allgatherv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize * getCommSize(h->comm));

    AMPI_Allgatherv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                            h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->comm);

    adjointInterface->combineAdjoints(h->sendbufAdjoints, h->sendbufTotalSize, getCommSize(h->comm));
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Allgatherv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                      typename RECVTYPE::Type* recvbuf, MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* displs,
                      RECVTYPE* recvtype, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Allgatherv(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcounts, displs,
                               recvtype->getMpiType(), comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* displsMod = displs;
      int displsTotalSize = 0;
      if(nullptr != displs) {
        displsTotalSize = computeDisplacementsTotalSize(recvcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          displsMod = createLinearDisplacements(recvcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcounts[getCommRank(comm)];
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = displsTotalSize;

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Allgatherv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->copyIntoModifiedBuffer(recvbuf, displs[rank], recvbufMod, displsMod[rank], recvcounts[rank]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(displs[getCommRank(comm)] + recvcounts[getCommRank(
                              comm)]) - recvtype->computeActiveElements(displs[getCommRank(comm)]);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        createLinearIndexCounts(h->recvbufCount, recvcounts, displs, getCommSize(comm), recvtype);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->getValues(recvbuf, displs[i], h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->getIndices(recvbuf, displs[rank], h->sendbufIndices, 0, recvcounts[rank]);
          }
        }

        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->createIndices(recvbuf, displs[i], h->recvbufIndices, displsMod[i], recvcounts[i]);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Allgatherv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Allgatherv_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Allgatherv_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcounts = recvcounts;
        h->displs = displs;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->clearIndices(recvbuf, displs[i], recvcounts[i]);
        }
      }

      rStatus = MPI_Allgatherv(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcounts, displsMod,
                               recvtype->getModifiedMpiType(), comm);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->copyFromModifiedBuffer(recvbuf, displs[i], recvbufMod, displsMod[i], recvcounts[i]);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->registerValue(recvbuf, displs[i], h->recvbufIndices, h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] displsMod;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Allreduce_global_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename DATATYPE::IndexType* sendbufIndices;
    typename DATATYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int recvbufTotalSize;
    typename DATATYPE::IndexType* recvbufIndices;
    typename DATATYPE::PrimalType* recvbufPrimals;
    typename DATATYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int count;
    DATATYPE* datatype;
    AMPI_Op op;
    AMPI_Comm comm;

    ~AMPI_Allreduce_global_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Allreduce_global_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Allreduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Allreduce_global_pri<DATATYPE>(h->sendbufPrimals, h->sendbufCountVec, h->recvbufPrimals, h->recvbufCountVec,
                                        h->count, h->datatype, h->op, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Allreduce_global_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Allreduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Allreduce_global_fwd<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                        h->count, h->datatype, h->op, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Allreduce_global_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Allreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Allreduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    convOp.preAdjointOperation(h->recvbufAdjoints, h->recvbufPrimals, h->recvbufCount, adjointInterface->getVectorSize());
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize * getCommSize(h->comm));

    AMPI_Allreduce_global_adj<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                        h->count, h->datatype, h->op, h->comm);

    adjointInterface->combineAdjoints(h->sendbufAdjoints, h->sendbufTotalSize, getCommSize(h->comm));
    // the primals of the recive buffer are always given to the function. The operator should ignore them if not needed.
    // The wrapper functions make sure that for operators that need the primals an all* action is perfomed (e.g. Allreduce instead of Reduce)
    convOp.postAdjointOperation(h->sendbufAdjoints, h->sendbufPrimals, h->recvbufPrimals, h->sendbufTotalSize,
                                adjointInterface->getVectorSize());
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Allreduce_global(MEDI_OPTIONAL_CONST typename DATATYPE::Type* sendbuf, typename DATATYPE::Type* recvbuf,
                            int count, DATATYPE* datatype, AMPI_Op op, AMPI_Comm comm) {
    int rStatus;
    AMPI_Op convOp = datatype->getADTool().convertOperator(op);
    (void)convOp;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Allreduce(sendbuf, recvbuf, count, datatype->getMpiType(), convOp.primalFunction, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Allreduce_global_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = count;
      } else {
        sendbufElements = count;
      }

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(sendbuf));
      }
      typename DATATYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Allreduce_global_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          datatype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, count);
        } else {
          datatype->copyIntoModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = datatype->computeActiveElements(count);
        } else {
          h->sendbufCount = datatype->computeActiveElements(count);
        }
        h->sendbufTotalSize = datatype->computeActiveElements(sendbufElements);
        datatype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = datatype->computeActiveElements(count);
        h->recvbufTotalSize = datatype->computeActiveElements(recvbufElements);
        datatype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);

        // extract the primal values for the operator if required
        if(convOp.requiresPrimal) {
          datatype->getADTool().createPrimalTypeBuffer(h->sendbufPrimals, h->sendbufTotalSize);
          if(AMPI_IN_PLACE != sendbuf) {
            datatype->getValues(sendbuf, 0, h->sendbufPrimals, 0, count);
          } else {
            datatype->getValues(recvbuf, 0, h->sendbufPrimals, 0, count);
          }
        }

        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          datatype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, count);
        }


        if(AMPI_IN_PLACE != sendbuf) {
          datatype->getIndices(sendbuf, 0, h->sendbufIndices, 0, count);
        } else {
          datatype->getIndices(recvbuf, 0, h->sendbufIndices, 0, count);
        }

        datatype->createIndices(recvbuf, 0, h->recvbufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Allreduce_global_b<DATATYPE>;
        h->funcForward = AMPI_Allreduce_global_d<DATATYPE>;
        h->funcPrimal = AMPI_Allreduce_global_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->op = op;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(recvbuf, 0, count);
      }

      rStatus = MPI_Allreduce(sendbufMod, recvbufMod, count, datatype->getModifiedMpiType(), convOp.modifiedPrimalFunction,
                              comm);
      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, count);
      }
      // extract the primal values for the operator if required
      if(nullptr != h && convOp.requiresPrimal) {
        datatype->getADTool().createPrimalTypeBuffer(h->recvbufPrimals, h->recvbufTotalSize);
        datatype->getValues(recvbuf, 0, h->recvbufPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Alltoall_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Alltoall_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoall_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Alltoall_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoall_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Alltoall_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoall_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Alltoall_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Alltoall(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                    typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Alltoall(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount * getCommSize(comm);
      } else {
        sendbufElements = recvcount * getCommSize(comm);
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = recvcount * getCommSize(comm);

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Alltoall_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount * getCommSize(comm));
        } else {
          recvtype->copyIntoModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(recvcount);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = recvtype->computeActiveElements(recvcount);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount * getCommSize(comm));
        } else {
          recvtype->getIndices(recvbuf, 0, h->sendbufIndices, 0, recvcount * getCommSize(comm));
        }

        recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount * getCommSize(comm));

        // pack all the variables in the handle
        h->funcReverse = AMPI_Alltoall_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Alltoall_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Alltoall_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        recvtype->clearIndices(recvbuf, 0, recvcount * getCommSize(comm));
      }

      rStatus = MPI_Alltoall(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                             recvtype->getModifiedMpiType(), comm);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
      }

      if(nullptr != h) {
        // handle the recv buffers
        recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Alltoallv_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int* sendbufCount;
    /* required for async */ int* sendbufCountVec;
    /* required for async */ int* sendbufDisplsVec;
    MEDI_OPTIONAL_CONST  int* sendcounts;
    MEDI_OPTIONAL_CONST  int* sdispls;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int* recvbufCount;
    /* required for async */ int* recvbufCountVec;
    /* required for async */ int* recvbufDisplsVec;
    MEDI_OPTIONAL_CONST  int* recvcounts;
    MEDI_OPTIONAL_CONST  int* rdispls;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Alltoallv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != sendbufCount) {
        delete [] sendbufCount;
        sendbufCount = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
      if(nullptr != recvbufCount) {
        delete [] recvbufCount;
        recvbufCount = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoallv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Alltoallv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                           h->sdispls, h->sendtype, h->recvbufPrimals, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->rdispls,
                                           h->recvtype, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    delete [] h->sendbufCountVec;
    delete [] h->sendbufDisplsVec;
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoallv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Alltoallv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                           h->sdispls, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->rdispls,
                                           h->recvtype, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    delete [] h->sendbufCountVec;
    delete [] h->sendbufDisplsVec;
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Alltoallv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Alltoallv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                           h->sdispls, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->rdispls,
                                           h->recvtype, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    delete [] h->sendbufCountVec;
    delete [] h->sendbufDisplsVec;
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Alltoallv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, MEDI_OPTIONAL_CONST int* sendcounts,
                     MEDI_OPTIONAL_CONST int* sdispls, SENDTYPE* sendtype, typename RECVTYPE::Type* recvbuf,
                     MEDI_OPTIONAL_CONST int* recvcounts, MEDI_OPTIONAL_CONST int* rdispls, RECVTYPE* recvtype, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype->getMpiType(), recvbuf, recvcounts, rdispls,
                              recvtype->getMpiType(), comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* sdisplsMod = sdispls;
      int sdisplsTotalSize = 0;
      if(nullptr != sdispls) {
        sdisplsTotalSize = computeDisplacementsTotalSize(sendcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          sdisplsMod = createLinearDisplacements(sendcounts, getCommSize(comm));
        }
      }
      MEDI_OPTIONAL_CONST int* rdisplsMod = rdispls;
      int rdisplsTotalSize = 0;
      if(nullptr != rdispls) {
        rdisplsTotalSize = computeDisplacementsTotalSize(recvcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          rdisplsMod = createLinearDisplacements(recvcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sdisplsTotalSize;
      } else {
        sendbufElements = rdisplsTotalSize;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = rdisplsTotalSize;

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Alltoallv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->copyIntoModifiedBuffer(sendbuf, sdispls[i], sendbufMod, sdisplsMod[i], sendcounts[i]);
          }
        } else {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->copyIntoModifiedBuffer(recvbuf, rdispls[i], recvbufMod, rdisplsMod[i], recvcounts[i]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          createLinearIndexCounts(h->sendbufCount, sendcounts, sdispls, getCommSize(comm), sendtype);
        } else {
          createLinearIndexCounts(h->sendbufCount, recvcounts, rdispls, getCommSize(comm), recvtype);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        createLinearIndexCounts(h->recvbufCount, recvcounts, rdispls, getCommSize(comm), recvtype);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->getValues(recvbuf, rdispls[i], h->recvbufOldPrimals, rdisplsMod[i], recvcounts[i]);
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->getIndices(sendbuf, sdispls[i], h->sendbufIndices, sdisplsMod[i], sendcounts[i]);
          }
        } else {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->getIndices(recvbuf, rdispls[i], h->sendbufIndices, rdisplsMod[i], recvcounts[i]);
          }
        }

        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->createIndices(recvbuf, rdispls[i], h->recvbufIndices, rdisplsMod[i], recvcounts[i]);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Alltoallv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Alltoallv_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Alltoallv_p<SENDTYPE, RECVTYPE>;
        h->sendcounts = sendcounts;
        h->sdispls = sdispls;
        h->sendtype = sendtype;
        h->recvcounts = recvcounts;
        h->rdispls = rdispls;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->clearIndices(recvbuf, rdispls[i], recvcounts[i]);
        }
      }

      rStatus = MPI_Alltoallv(sendbufMod, sendcounts, sdisplsMod, sendtype->getModifiedMpiType(), recvbufMod, recvcounts,
                              rdisplsMod, recvtype->getModifiedMpiType(), comm);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->copyFromModifiedBuffer(recvbuf, rdispls[i], recvbufMod, rdisplsMod[i], recvcounts[i]);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->registerValue(recvbuf, rdispls[i], h->recvbufIndices, h->recvbufOldPrimals, rdisplsMod[i], recvcounts[i]);
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] sdisplsMod;
      }
      if(recvtype->isModifiedBufferRequired()) {
        delete [] rdisplsMod;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Bcast_wrap_AdjointHandle : public HandleBase {
    int bufferSendTotalSize;
    typename DATATYPE::IndexType* bufferSendIndices;
    typename DATATYPE::PrimalType* bufferSendPrimals;
    /* required for async */ void* bufferSendAdjoints;
    int bufferSendCount;
    int bufferSendCountVec;
    int bufferRecvTotalSize;
    typename DATATYPE::IndexType* bufferRecvIndices;
    typename DATATYPE::PrimalType* bufferRecvPrimals;
    typename DATATYPE::PrimalType* bufferRecvOldPrimals;
    /* required for async */ void* bufferRecvAdjoints;
    int bufferRecvCount;
    int bufferRecvCountVec;
    int count;
    DATATYPE* datatype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Bcast_wrap_AdjointHandle () {
      if(nullptr != bufferSendIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufferSendIndices);
        bufferSendIndices = nullptr;
      }
      if(nullptr != bufferSendPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufferSendPrimals);
        bufferSendPrimals = nullptr;
      }
      if(nullptr != bufferRecvIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufferRecvIndices);
        bufferRecvIndices = nullptr;
      }
      if(nullptr != bufferRecvPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufferRecvPrimals);
        bufferRecvPrimals = nullptr;
      }
      if(nullptr != bufferRecvOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufferRecvOldPrimals);
        bufferRecvOldPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Bcast_wrap_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Bcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Bcast_wrap_AdjointHandle<DATATYPE>*>(handle);

    h->bufferRecvAdjoints = nullptr;
    h->bufferRecvCountVec = adjointInterface->getVectorSize() * h->bufferRecvCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufferRecvPrimals, h->bufferRecvTotalSize );
    h->bufferSendAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->bufferSendCountVec = adjointInterface->getVectorSize() * h->bufferSendCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->bufferSendPrimals, h->bufferSendTotalSize );
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getPrimals(h->bufferSendIndices, h->bufferSendPrimals, h->bufferSendTotalSize);

    }

    AMPI_Bcast_wrap_pri<DATATYPE>(h->bufferSendPrimals, h->bufferSendCountVec, h->bufferRecvPrimals, h->bufferRecvCountVec,
                                  h->count, h->datatype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deletePrimalTypeBuffer((void*&)h->bufferSendPrimals);
    }
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->bufferRecvIndices, h->bufferRecvOldPrimals, h->bufferRecvTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->bufferRecvIndices, h->bufferRecvPrimals, h->bufferRecvTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufferRecvPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Bcast_wrap_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Bcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Bcast_wrap_AdjointHandle<DATATYPE>*>(handle);

    h->bufferRecvAdjoints = nullptr;
    h->bufferRecvCountVec = adjointInterface->getVectorSize() * h->bufferRecvCount;
    adjointInterface->createAdjointTypeBuffer(h->bufferRecvAdjoints, h->bufferRecvTotalSize );
    h->bufferSendAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->bufferSendCountVec = adjointInterface->getVectorSize() * h->bufferSendCount;
      adjointInterface->createAdjointTypeBuffer(h->bufferSendAdjoints, h->bufferSendTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->bufferSendIndices, h->bufferSendAdjoints, h->bufferSendTotalSize);

    }

    AMPI_Bcast_wrap_fwd<DATATYPE>(h->bufferSendAdjoints, h->bufferSendCountVec, h->bufferRecvAdjoints,
                                  h->bufferRecvCountVec, h->count, h->datatype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->bufferSendAdjoints);
    }
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufferRecvIndices, h->bufferRecvAdjoints, h->bufferRecvTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufferRecvAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Bcast_wrap_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Bcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Bcast_wrap_AdjointHandle<DATATYPE>*>(handle);

    h->bufferRecvAdjoints = nullptr;
    h->bufferRecvCountVec = adjointInterface->getVectorSize() * h->bufferRecvCount;
    adjointInterface->createAdjointTypeBuffer(h->bufferRecvAdjoints, h->bufferRecvTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufferRecvIndices, h->bufferRecvAdjoints, h->bufferRecvTotalSize);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->bufferRecvIndices, h->bufferRecvOldPrimals, h->bufferRecvTotalSize);
    }
    h->bufferSendAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->bufferSendCountVec = adjointInterface->getVectorSize() * h->bufferSendCount;
      adjointInterface->createAdjointTypeBuffer(h->bufferSendAdjoints, h->bufferSendTotalSize * getCommSize(h->comm));
    }

    AMPI_Bcast_wrap_adj<DATATYPE>(h->bufferSendAdjoints, h->bufferSendCountVec, h->bufferRecvAdjoints,
                                  h->bufferRecvCountVec, h->count, h->datatype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->combineAdjoints(h->bufferSendAdjoints, h->bufferSendTotalSize, getCommSize(h->comm));
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->bufferSendIndices, h->bufferSendAdjoints, h->bufferSendTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->bufferSendAdjoints);
    }
    adjointInterface->deleteAdjointTypeBuffer(h->bufferRecvAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Bcast_wrap(typename DATATYPE::Type* bufferSend, typename DATATYPE::Type* bufferRecv, int count,
                      DATATYPE* datatype, int root, AMPI_Comm comm) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Bcast_wrap(bufferSend, bufferRecv, count, datatype->getMpiType(), root, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Bcast_wrap_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufferSendMod = nullptr;
      int bufferSendElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        if(AMPI_IN_PLACE != bufferSend) {
          bufferSendElements = count;
        } else {
          bufferSendElements = count;
        }

        if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == bufferSend)) {
          datatype->createModifiedTypeBuffer(bufferSendMod, bufferSendElements);
        } else {
          bufferSendMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(bufferSend));
        }
      }
      typename DATATYPE::ModifiedType* bufferRecvMod = nullptr;
      int bufferRecvElements = 0;

      // compute the total size of the buffer
      bufferRecvElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufferRecvMod, bufferRecvElements);
      } else {
        bufferRecvMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(bufferRecv));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Bcast_wrap_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired()) {
          if(AMPI_IN_PLACE != bufferSend) {
            datatype->copyIntoModifiedBuffer(bufferSend, 0, bufferSendMod, 0, count);
          } else {
            datatype->copyIntoModifiedBuffer(bufferRecv, 0, bufferRecvMod, 0, count);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(root == getCommRank(comm)) {
          if(AMPI_IN_PLACE != bufferSend) {
            h->bufferSendCount = datatype->computeActiveElements(count);
          } else {
            h->bufferSendCount = datatype->computeActiveElements(count);
          }
          h->bufferSendTotalSize = datatype->computeActiveElements(bufferSendElements);
          datatype->getADTool().createIndexTypeBuffer(h->bufferSendIndices, h->bufferSendTotalSize);
        }
        h->bufferRecvCount = datatype->computeActiveElements(count);
        h->bufferRecvTotalSize = datatype->computeActiveElements(bufferRecvElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufferRecvIndices, h->bufferRecvTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufferRecvOldPrimals, h->bufferRecvTotalSize);
          datatype->getValues(bufferRecv, 0, h->bufferRecvOldPrimals, 0, count);
        }


        if(root == getCommRank(comm)) {
          if(AMPI_IN_PLACE != bufferSend) {
            datatype->getIndices(bufferSend, 0, h->bufferSendIndices, 0, count);
          } else {
            datatype->getIndices(bufferRecv, 0, h->bufferSendIndices, 0, count);
          }
        }

        datatype->createIndices(bufferRecv, 0, h->bufferRecvIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Bcast_wrap_b<DATATYPE>;
        h->funcForward = AMPI_Bcast_wrap_d<DATATYPE>;
        h->funcPrimal = AMPI_Bcast_wrap_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->root = root;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(bufferRecv, 0, count);
      }

      rStatus = MPI_Bcast_wrap(bufferSendMod, bufferRecvMod, count, datatype->getModifiedMpiType(), root, comm);
      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(bufferRecv, 0, bufferRecvMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(bufferRecv, 0, h->bufferRecvIndices, h->bufferRecvOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == bufferSend)) {
          datatype->deleteModifiedTypeBuffer(bufferSendMod);
        }
      }
      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufferRecvMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Gather_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Gather_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gather_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Gather_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype, h->recvbufPrimals,
                                        h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    if(h->root == getCommRank(h->comm)) {
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
      adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gather_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Gather_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                        h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gather_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Gather_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                        h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Gather(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                  typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Gather(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), root,
                           comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcount;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        recvbufElements = recvcount * getCommSize(comm);

        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
        } else {
          recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
        }
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Gather_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          recvtype->copyIntoModifiedBuffer(recvbuf, recvcount * getCommRank(comm), recvbufMod, recvcount * getCommRank(comm),
                                           recvcount);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(recvcount);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        if(root == getCommRank(comm)) {
          h->recvbufCount = recvtype->computeActiveElements(recvcount);
          h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);
        }


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          if(root == getCommRank(comm)) {
            recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
            if(root == getCommRank(comm)) {
              recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
            }
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          recvtype->getIndices(recvbuf, recvcount * getCommRank(comm), h->sendbufIndices, 0, recvcount);
        }

        if(root == getCommRank(comm)) {
          recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount * getCommSize(comm));
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Gather_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Gather_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Gather_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(root == getCommRank(comm)) {
        if(!recvtype->isModifiedBufferRequired()) {
          recvtype->clearIndices(recvbuf, 0, recvcount * getCommSize(comm));
        }
      }

      rStatus = MPI_Gather(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                           recvtype->getModifiedMpiType(), root, comm);
      recvtype->getADTool().addToolAction(h);

      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired()) {
          recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(root == getCommRank(comm)) {
          recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
        }
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->deleteModifiedTypeBuffer(recvbufMod);
        }
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Gatherv_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int* recvbufCount;
    /* required for async */ int* recvbufCountVec;
    /* required for async */ int* recvbufDisplsVec;
    const  int* recvcounts;
    const  int* displs;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Gatherv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
      if(nullptr != recvbufCount) {
        delete [] recvbufCount;
        recvbufCount = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gatherv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Gatherv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufPrimals, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->root, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    if(h->root == getCommRank(h->comm)) {
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
      adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
      delete [] h->recvbufCountVec;
      delete [] h->recvbufDisplsVec;
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gatherv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Gatherv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->root, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
      delete [] h->recvbufCountVec;
      delete [] h->recvbufDisplsVec;
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Gatherv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Gatherv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->root, h->comm);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
      delete [] h->recvbufCountVec;
      delete [] h->recvbufDisplsVec;
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Gatherv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                   typename RECVTYPE::Type* recvbuf, const int* recvcounts, const int* displs, RECVTYPE* recvtype, int root,
                   AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Gatherv(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcounts, displs, recvtype->getMpiType(),
                            root, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* displsMod = displs;
      int displsTotalSize = 0;
      if(nullptr != displs) {
        displsTotalSize = computeDisplacementsTotalSize(recvcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          displsMod = createLinearDisplacements(recvcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcounts[getCommRank(comm)];
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        recvbufElements = displsTotalSize;

        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
        } else {
          recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
        }
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Gatherv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->copyIntoModifiedBuffer(recvbuf, displs[rank], recvbufMod, displsMod[rank], recvcounts[rank]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(displs[getCommRank(comm)] + recvcounts[getCommRank(
                              comm)]) - recvtype->computeActiveElements(displs[getCommRank(comm)]);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        if(root == getCommRank(comm)) {
          createLinearIndexCounts(h->recvbufCount, recvcounts, displs, getCommSize(comm), recvtype);
          h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);
        }


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          if(root == getCommRank(comm)) {
            recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
            if(root == getCommRank(comm)) {
              for(int i = 0; i < getCommSize(comm); ++i) {
                recvtype->getValues(recvbuf, displs[i], h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
              }
            }
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->getIndices(recvbuf, displs[rank], h->sendbufIndices, 0, recvcounts[rank]);
          }
        }

        if(root == getCommRank(comm)) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->createIndices(recvbuf, displs[i], h->recvbufIndices, displsMod[i], recvcounts[i]);
          }
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Gatherv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Gatherv_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Gatherv_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcounts = recvcounts;
        h->displs = displs;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(root == getCommRank(comm)) {
        if(!recvtype->isModifiedBufferRequired()) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->clearIndices(recvbuf, displs[i], recvcounts[i]);
          }
        }
      }

      rStatus = MPI_Gatherv(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcounts, displsMod,
                            recvtype->getModifiedMpiType(), root, comm);
      recvtype->getADTool().addToolAction(h);

      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired()) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->copyFromModifiedBuffer(recvbuf, displs[i], recvbufMod, displsMod[i], recvcounts[i]);
          }
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(root == getCommRank(comm)) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->registerValue(recvbuf, displs[i], h->recvbufIndices, h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
          }
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] displsMod;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->deleteModifiedTypeBuffer(recvbufMod);
        }
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iallgather_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Iallgather_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iallgather_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    int sendcount;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    int recvcount;
    RECVTYPE* recvtype;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Iallgather_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                            h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Iallgather_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                            h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize * getCommSize(h->comm));

    AMPI_Iallgather_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                            h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgather_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->combineAdjoints(h->sendbufAdjoints, h->sendbufTotalSize, getCommSize(h->comm));
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iallgather_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iallgather(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                      typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Iallgather(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), comm,
                               &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcount;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = recvcount * getCommSize(comm);

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          recvtype->copyIntoModifiedBuffer(recvbuf, recvcount * getCommRank(comm), recvbufMod, recvcount * getCommRank(comm),
                                           recvcount);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(recvcount);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = recvtype->computeActiveElements(recvcount);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          recvtype->getIndices(recvbuf, recvcount * getCommRank(comm), h->sendbufIndices, 0, recvcount);
        }

        recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount * getCommSize(comm));

        // pack all the variables in the handle
        h->funcReverse = AMPI_Iallgather_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Iallgather_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Iallgather_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        recvtype->clearIndices(recvbuf, 0, recvcount * getCommSize(comm));
      }

      rStatus = MPI_Iallgather(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                               recvtype->getModifiedMpiType(), comm, &request->request);

      AMPI_Iallgather_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Iallgather_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sendcount = sendcount;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->recvcount = recvcount;
      asyncHandle->recvtype = recvtype;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Iallgather_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Iallgather_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Iallgather_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iallgather_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Iallgather_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Iallgather_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    int sendcount = asyncHandle->sendcount;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int recvcount = asyncHandle->recvcount;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iallgather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcount); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcount); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
      }

      if(nullptr != h) {
        // handle the recv buffers
        recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iallgatherv_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int* recvbufCount;
    /* required for async */ int* recvbufCountVec;
    /* required for async */ int* recvbufDisplsVec;
    const  int* recvcounts;
    const  int* displs;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Iallgatherv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
      if(nullptr != recvbufCount) {
        delete [] recvbufCount;
        recvbufCount = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iallgatherv_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    int sendcount;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    const int* displsMod;
    const  int* recvcounts;
    const  int* displs;
    RECVTYPE* recvtype;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Iallgatherv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
        h->recvbufPrimals, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->comm,
        &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Iallgatherv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
        h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->comm,
        &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize * getCommSize(h->comm));

    AMPI_Iallgatherv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
        h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->comm,
        &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iallgatherv_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->combineAdjoints(h->sendbufAdjoints, h->sendbufTotalSize, getCommSize(h->comm));
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iallgatherv_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iallgatherv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                       typename RECVTYPE::Type* recvbuf, const int* recvcounts, const int* displs, RECVTYPE* recvtype, AMPI_Comm comm,
                       AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Iallgatherv(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcounts, displs,
                                recvtype->getMpiType(), comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* displsMod = displs;
      int displsTotalSize = 0;
      if(nullptr != displs) {
        displsTotalSize = computeDisplacementsTotalSize(recvcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          displsMod = createLinearDisplacements(recvcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcounts[getCommRank(comm)];
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = displsTotalSize;

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->copyIntoModifiedBuffer(recvbuf, displs[rank], recvbufMod, displsMod[rank], recvcounts[rank]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(displs[getCommRank(comm)] + recvcounts[getCommRank(
                              comm)]) - recvtype->computeActiveElements(displs[getCommRank(comm)]);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        createLinearIndexCounts(h->recvbufCount, recvcounts, displs, getCommSize(comm), recvtype);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->getValues(recvbuf, displs[i], h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->getIndices(recvbuf, displs[rank], h->sendbufIndices, 0, recvcounts[rank]);
          }
        }

        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->createIndices(recvbuf, displs[i], h->recvbufIndices, displsMod[i], recvcounts[i]);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Iallgatherv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Iallgatherv_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Iallgatherv_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcounts = recvcounts;
        h->displs = displs;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->clearIndices(recvbuf, displs[i], recvcounts[i]);
        }
      }

      rStatus = MPI_Iallgatherv(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcounts, displsMod,
                                recvtype->getModifiedMpiType(), comm, &request->request);

      AMPI_Iallgatherv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Iallgatherv_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sendcount = sendcount;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->displsMod = displsMod;
      asyncHandle->recvcounts = recvcounts;
      asyncHandle->displs = displs;
      asyncHandle->recvtype = recvtype;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Iallgatherv_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Iallgatherv_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Iallgatherv_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iallgatherv_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Iallgatherv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Iallgatherv_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    int sendcount = asyncHandle->sendcount;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    const int* displsMod = asyncHandle->displsMod;
    const  int* recvcounts = asyncHandle->recvcounts;
    const  int* displs = asyncHandle->displs;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h =
      static_cast<AMPI_Iallgatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>(asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcount); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(displsMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcounts); // Unused generated to ignore warnings
    MEDI_UNUSED(displs); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->copyFromModifiedBuffer(recvbuf, displs[i], recvbufMod, displsMod[i], recvcounts[i]);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->registerValue(recvbuf, displs[i], h->recvbufIndices, h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] displsMod;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Iallreduce_global_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename DATATYPE::IndexType* sendbufIndices;
    typename DATATYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int recvbufTotalSize;
    typename DATATYPE::IndexType* recvbufIndices;
    typename DATATYPE::PrimalType* recvbufPrimals;
    typename DATATYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int count;
    DATATYPE* datatype;
    AMPI_Op op;
    AMPI_Comm comm;

    ~AMPI_Iallreduce_global_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Iallreduce_global_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* sendbuf;
    typename DATATYPE::ModifiedType* sendbufMod;
    typename DATATYPE::Type* recvbuf;
    typename DATATYPE::ModifiedType* recvbufMod;
    int count;
    DATATYPE* datatype;
    AMPI_Op op;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Iallreduce_global_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Iallreduce_global_pri<DATATYPE>(h->sendbufPrimals, h->sendbufCountVec, h->recvbufPrimals, h->recvbufCountVec,
                                         h->count, h->datatype, h->op, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Iallreduce_global_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Iallreduce_global_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Iallreduce_global_fwd<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                         h->count, h->datatype, h->op, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Iallreduce_global_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Iallreduce_global_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    convOp.preAdjointOperation(h->recvbufAdjoints, h->recvbufPrimals, h->recvbufCount, adjointInterface->getVectorSize());
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize * getCommSize(h->comm));

    AMPI_Iallreduce_global_adj<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                         h->count, h->datatype, h->op, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Iallreduce_global_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    adjointInterface->combineAdjoints(h->sendbufAdjoints, h->sendbufTotalSize, getCommSize(h->comm));
    // the primals of the recive buffer are always given to the function. The operator should ignore them if not needed.
    // The wrapper functions make sure that for operators that need the primals an all* action is perfomed (e.g. Allreduce instead of Reduce)
    convOp.postAdjointOperation(h->sendbufAdjoints, h->sendbufPrimals, h->recvbufPrimals, h->sendbufTotalSize,
                                adjointInterface->getVectorSize());
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Iallreduce_global_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Iallreduce_global(MEDI_OPTIONAL_CONST typename DATATYPE::Type* sendbuf, typename DATATYPE::Type* recvbuf,
                             int count, DATATYPE* datatype, AMPI_Op op, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;
    AMPI_Op convOp = datatype->getADTool().convertOperator(op);
    (void)convOp;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Iallreduce(sendbuf, recvbuf, count, datatype->getMpiType(), convOp.primalFunction, comm,
                               &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = count;
      } else {
        sendbufElements = count;
      }

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(sendbuf));
      }
      typename DATATYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Iallreduce_global_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          datatype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, count);
        } else {
          datatype->copyIntoModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = datatype->computeActiveElements(count);
        } else {
          h->sendbufCount = datatype->computeActiveElements(count);
        }
        h->sendbufTotalSize = datatype->computeActiveElements(sendbufElements);
        datatype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = datatype->computeActiveElements(count);
        h->recvbufTotalSize = datatype->computeActiveElements(recvbufElements);
        datatype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);

        // extract the primal values for the operator if required
        if(convOp.requiresPrimal) {
          datatype->getADTool().createPrimalTypeBuffer(h->sendbufPrimals, h->sendbufTotalSize);
          if(AMPI_IN_PLACE != sendbuf) {
            datatype->getValues(sendbuf, 0, h->sendbufPrimals, 0, count);
          } else {
            datatype->getValues(recvbuf, 0, h->sendbufPrimals, 0, count);
          }
        }

        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          datatype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, count);
        }


        if(AMPI_IN_PLACE != sendbuf) {
          datatype->getIndices(sendbuf, 0, h->sendbufIndices, 0, count);
        } else {
          datatype->getIndices(recvbuf, 0, h->sendbufIndices, 0, count);
        }

        datatype->createIndices(recvbuf, 0, h->recvbufIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Iallreduce_global_b<DATATYPE>;
        h->funcForward = AMPI_Iallreduce_global_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Iallreduce_global_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->op = op;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(recvbuf, 0, count);
      }

      rStatus = MPI_Iallreduce(sendbufMod, recvbufMod, count, datatype->getModifiedMpiType(), convOp.modifiedPrimalFunction,
                               comm, &request->request);

      AMPI_Iallreduce_global_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Iallreduce_global_AsyncHandle<DATATYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->op = op;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Iallreduce_global_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Iallreduce_global_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Iallreduce_global_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Iallreduce_global_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Iallreduce_global_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Iallreduce_global_AsyncHandle<DATATYPE>*>
        (handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename DATATYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    typename DATATYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename DATATYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    AMPI_Op op = asyncHandle->op;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Iallreduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Iallreduce_global_AdjointHandle<DATATYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(op); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      AMPI_Op convOp = datatype->getADTool().convertOperator(op);
      (void)convOp;
      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, count);
      }
      // extract the primal values for the operator if required
      if(nullptr != h && convOp.requiresPrimal) {
        datatype->getADTool().createPrimalTypeBuffer(h->recvbufPrimals, h->recvbufTotalSize);
        datatype->getValues(recvbuf, 0, h->recvbufPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Ialltoall_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Ialltoall_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Ialltoall_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    int sendcount;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    int recvcount;
    RECVTYPE* recvtype;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Ialltoall_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                           h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Ialltoall_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                           h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Ialltoall_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                           h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoall_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Ialltoall_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Ialltoall(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                     typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ialltoall(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), comm,
                              &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount * getCommSize(comm);
      } else {
        sendbufElements = recvcount * getCommSize(comm);
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = recvcount * getCommSize(comm);

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount * getCommSize(comm));
        } else {
          recvtype->copyIntoModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(recvcount);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        h->recvbufCount = recvtype->computeActiveElements(recvcount);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount * getCommSize(comm));
        } else {
          recvtype->getIndices(recvbuf, 0, h->sendbufIndices, 0, recvcount * getCommSize(comm));
        }

        recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount * getCommSize(comm));

        // pack all the variables in the handle
        h->funcReverse = AMPI_Ialltoall_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Ialltoall_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Ialltoall_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        recvtype->clearIndices(recvbuf, 0, recvcount * getCommSize(comm));
      }

      rStatus = MPI_Ialltoall(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                              recvtype->getModifiedMpiType(), comm, &request->request);

      AMPI_Ialltoall_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Ialltoall_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sendcount = sendcount;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->recvcount = recvcount;
      asyncHandle->recvtype = recvtype;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Ialltoall_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Ialltoall_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Ialltoall_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Ialltoall_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ialltoall_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Ialltoall_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    int sendcount = asyncHandle->sendcount;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int recvcount = asyncHandle->recvcount;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoall_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcount); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcount); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
      }

      if(nullptr != h) {
        // handle the recv buffers
        recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Ialltoallv_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int* sendbufCount;
    /* required for async */ int* sendbufCountVec;
    /* required for async */ int* sendbufDisplsVec;
    const  int* sendcounts;
    const  int* sdispls;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int* recvbufCount;
    /* required for async */ int* recvbufCountVec;
    /* required for async */ int* recvbufDisplsVec;
    const  int* recvcounts;
    const  int* rdispls;
    RECVTYPE* recvtype;
    AMPI_Comm comm;

    ~AMPI_Ialltoallv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != sendbufCount) {
        delete [] sendbufCount;
        sendbufCount = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
      if(nullptr != recvbufCount) {
        delete [] recvbufCount;
        recvbufCount = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Ialltoallv_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    const int* sdisplsMod;
    const  int* sendcounts;
    const  int* sdispls;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    const int* rdisplsMod;
    const  int* recvcounts;
    const  int* rdispls;
    RECVTYPE* recvtype;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Ialltoallv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                            h->sdispls, h->sendtype, h->recvbufPrimals, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->rdispls,
                                            h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    delete [] h->sendbufCountVec;
    delete [] h->sendbufDisplsVec;
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Ialltoallv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                            h->sdispls, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->rdispls,
                                            h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    delete [] h->sendbufCountVec;
    delete [] h->sendbufDisplsVec;
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                      adjointInterface->getVectorSize());
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Ialltoallv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                            h->sdispls, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->rdispls,
                                            h->recvtype, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Ialltoallv_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    delete [] h->sendbufCountVec;
    delete [] h->sendbufDisplsVec;
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    delete [] h->recvbufCountVec;
    delete [] h->recvbufDisplsVec;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Ialltoallv_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Ialltoallv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, const int* sendcounts, const int* sdispls,
                      SENDTYPE* sendtype, typename RECVTYPE::Type* recvbuf, const int* recvcounts, const int* rdispls, RECVTYPE* recvtype,
                      AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype->getMpiType(), recvbuf, recvcounts, rdispls,
                               recvtype->getMpiType(), comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* sdisplsMod = sdispls;
      int sdisplsTotalSize = 0;
      if(nullptr != sdispls) {
        sdisplsTotalSize = computeDisplacementsTotalSize(sendcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          sdisplsMod = createLinearDisplacements(sendcounts, getCommSize(comm));
        }
      }
      MEDI_OPTIONAL_CONST int* rdisplsMod = rdispls;
      int rdisplsTotalSize = 0;
      if(nullptr != rdispls) {
        rdisplsTotalSize = computeDisplacementsTotalSize(recvcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          rdisplsMod = createLinearDisplacements(recvcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sdisplsTotalSize;
      } else {
        sendbufElements = rdisplsTotalSize;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      recvbufElements = rdisplsTotalSize;

      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->copyIntoModifiedBuffer(sendbuf, sdispls[i], sendbufMod, sdisplsMod[i], sendcounts[i]);
          }
        } else {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->copyIntoModifiedBuffer(recvbuf, rdispls[i], recvbufMod, rdisplsMod[i], recvcounts[i]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          createLinearIndexCounts(h->sendbufCount, sendcounts, sdispls, getCommSize(comm), sendtype);
        } else {
          createLinearIndexCounts(h->sendbufCount, recvcounts, rdispls, getCommSize(comm), recvtype);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        createLinearIndexCounts(h->recvbufCount, recvcounts, rdispls, getCommSize(comm), recvtype);
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->getValues(recvbuf, rdispls[i], h->recvbufOldPrimals, rdisplsMod[i], recvcounts[i]);
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->getIndices(sendbuf, sdispls[i], h->sendbufIndices, sdisplsMod[i], sendcounts[i]);
          }
        } else {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->getIndices(recvbuf, rdispls[i], h->sendbufIndices, rdisplsMod[i], recvcounts[i]);
          }
        }

        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->createIndices(recvbuf, rdispls[i], h->recvbufIndices, rdisplsMod[i], recvcounts[i]);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Ialltoallv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Ialltoallv_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Ialltoallv_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcounts = sendcounts;
        h->sdispls = sdispls;
        h->sendtype = sendtype;
        h->recvcounts = recvcounts;
        h->rdispls = rdispls;
        h->recvtype = recvtype;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->clearIndices(recvbuf, rdispls[i], recvcounts[i]);
        }
      }

      rStatus = MPI_Ialltoallv(sendbufMod, sendcounts, sdisplsMod, sendtype->getModifiedMpiType(), recvbufMod, recvcounts,
                               rdisplsMod, recvtype->getModifiedMpiType(), comm, &request->request);

      AMPI_Ialltoallv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Ialltoallv_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sdisplsMod = sdisplsMod;
      asyncHandle->sendcounts = sendcounts;
      asyncHandle->sdispls = sdispls;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->rdisplsMod = rdisplsMod;
      asyncHandle->recvcounts = recvcounts;
      asyncHandle->rdispls = rdispls;
      asyncHandle->recvtype = recvtype;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Ialltoallv_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Ialltoallv_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Ialltoallv_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Ialltoallv_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ialltoallv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Ialltoallv_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    const int* sdisplsMod = asyncHandle->sdisplsMod;
    const  int* sendcounts = asyncHandle->sendcounts;
    const  int* sdispls = asyncHandle->sdispls;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    const int* rdisplsMod = asyncHandle->rdisplsMod;
    const  int* recvcounts = asyncHandle->recvcounts;
    const  int* rdispls = asyncHandle->rdispls;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Ialltoallv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sdisplsMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcounts); // Unused generated to ignore warnings
    MEDI_UNUSED(sdispls); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(rdisplsMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcounts); // Unused generated to ignore warnings
    MEDI_UNUSED(rdispls); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->copyFromModifiedBuffer(recvbuf, rdispls[i], recvbufMod, rdisplsMod[i], recvcounts[i]);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        for(int i = 0; i < getCommSize(comm); ++i) {
          recvtype->registerValue(recvbuf, rdispls[i], h->recvbufIndices, h->recvbufOldPrimals, rdisplsMod[i], recvcounts[i]);
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] sdisplsMod;
      }
      if(recvtype->isModifiedBufferRequired()) {
        delete [] rdisplsMod;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(recvtype->isModifiedBufferRequired() ) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Ibcast_wrap_AdjointHandle : public AsyncAdjointHandle {
    int bufferSendTotalSize;
    typename DATATYPE::IndexType* bufferSendIndices;
    typename DATATYPE::PrimalType* bufferSendPrimals;
    /* required for async */ void* bufferSendAdjoints;
    int bufferSendCount;
    int bufferSendCountVec;
    int bufferRecvTotalSize;
    typename DATATYPE::IndexType* bufferRecvIndices;
    typename DATATYPE::PrimalType* bufferRecvPrimals;
    typename DATATYPE::PrimalType* bufferRecvOldPrimals;
    /* required for async */ void* bufferRecvAdjoints;
    int bufferRecvCount;
    int bufferRecvCountVec;
    int count;
    DATATYPE* datatype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Ibcast_wrap_AdjointHandle () {
      if(nullptr != bufferSendIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufferSendIndices);
        bufferSendIndices = nullptr;
      }
      if(nullptr != bufferSendPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufferSendPrimals);
        bufferSendPrimals = nullptr;
      }
      if(nullptr != bufferRecvIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(bufferRecvIndices);
        bufferRecvIndices = nullptr;
      }
      if(nullptr != bufferRecvPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufferRecvPrimals);
        bufferRecvPrimals = nullptr;
      }
      if(nullptr != bufferRecvOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(bufferRecvOldPrimals);
        bufferRecvOldPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Ibcast_wrap_AsyncHandle : public AsyncHandle {
    typename DATATYPE::Type* bufferSend;
    typename DATATYPE::ModifiedType* bufferSendMod;
    typename DATATYPE::Type* bufferRecv;
    typename DATATYPE::ModifiedType* bufferRecvMod;
    int count;
    DATATYPE* datatype;
    int root;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>(handle);

    h->bufferRecvAdjoints = nullptr;
    h->bufferRecvCountVec = adjointInterface->getVectorSize() * h->bufferRecvCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->bufferRecvPrimals, h->bufferRecvTotalSize );
    h->bufferSendAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->bufferSendCountVec = adjointInterface->getVectorSize() * h->bufferSendCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->bufferSendPrimals, h->bufferSendTotalSize );
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getPrimals(h->bufferSendIndices, h->bufferSendPrimals, h->bufferSendTotalSize);

    }

    AMPI_Ibcast_wrap_pri<DATATYPE>(h->bufferSendPrimals, h->bufferSendCountVec, h->bufferRecvPrimals, h->bufferRecvCountVec,
                                   h->count, h->datatype, h->root, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deletePrimalTypeBuffer((void*&)h->bufferSendPrimals);
    }
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->bufferRecvIndices, h->bufferRecvOldPrimals, h->bufferRecvTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->bufferRecvIndices, h->bufferRecvPrimals, h->bufferRecvTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->bufferRecvPrimals);
  }

  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>(handle);

    h->bufferRecvAdjoints = nullptr;
    h->bufferRecvCountVec = adjointInterface->getVectorSize() * h->bufferRecvCount;
    adjointInterface->createAdjointTypeBuffer(h->bufferRecvAdjoints, h->bufferRecvTotalSize );
    h->bufferSendAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->bufferSendCountVec = adjointInterface->getVectorSize() * h->bufferSendCount;
      adjointInterface->createAdjointTypeBuffer(h->bufferSendAdjoints, h->bufferSendTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->bufferSendIndices, h->bufferSendAdjoints, h->bufferSendTotalSize);

    }

    AMPI_Ibcast_wrap_fwd<DATATYPE>(h->bufferSendAdjoints, h->bufferSendCountVec, h->bufferRecvAdjoints,
                                   h->bufferRecvCountVec, h->count, h->datatype, h->root, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->bufferSendAdjoints);
    }
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->bufferRecvIndices, h->bufferRecvAdjoints, h->bufferRecvTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->bufferRecvAdjoints);
  }

  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>(handle);

    h->bufferRecvAdjoints = nullptr;
    h->bufferRecvCountVec = adjointInterface->getVectorSize() * h->bufferRecvCount;
    adjointInterface->createAdjointTypeBuffer(h->bufferRecvAdjoints, h->bufferRecvTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->bufferRecvIndices, h->bufferRecvAdjoints, h->bufferRecvTotalSize);

    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->bufferRecvIndices, h->bufferRecvOldPrimals, h->bufferRecvTotalSize);
    }
    h->bufferSendAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->bufferSendCountVec = adjointInterface->getVectorSize() * h->bufferSendCount;
      adjointInterface->createAdjointTypeBuffer(h->bufferSendAdjoints, h->bufferSendTotalSize * getCommSize(h->comm));
    }

    AMPI_Ibcast_wrap_adj<DATATYPE>(h->bufferSendAdjoints, h->bufferSendCountVec, h->bufferRecvAdjoints,
                                   h->bufferRecvCountVec, h->count, h->datatype, h->root, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ibcast_wrap_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->combineAdjoints(h->bufferSendAdjoints, h->bufferSendTotalSize, getCommSize(h->comm));
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->bufferSendIndices, h->bufferSendAdjoints, h->bufferSendTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->bufferSendAdjoints);
    }
    adjointInterface->deleteAdjointTypeBuffer(h->bufferRecvAdjoints);
  }

  template<typename DATATYPE>
  int AMPI_Ibcast_wrap_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Ibcast_wrap(typename DATATYPE::Type* bufferSend, typename DATATYPE::Type* bufferRecv, int count,
                       DATATYPE* datatype, int root, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ibcast_wrap(bufferSend, bufferRecv, count, datatype->getMpiType(), root, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* bufferSendMod = nullptr;
      int bufferSendElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        if(AMPI_IN_PLACE != bufferSend) {
          bufferSendElements = count;
        } else {
          bufferSendElements = count;
        }

        if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == bufferSend)) {
          datatype->createModifiedTypeBuffer(bufferSendMod, bufferSendElements);
        } else {
          bufferSendMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(bufferSend));
        }
      }
      typename DATATYPE::ModifiedType* bufferRecvMod = nullptr;
      int bufferRecvElements = 0;

      // compute the total size of the buffer
      bufferRecvElements = count;

      if(datatype->isModifiedBufferRequired() ) {
        datatype->createModifiedTypeBuffer(bufferRecvMod, bufferRecvElements);
      } else {
        bufferRecvMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(bufferRecv));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired()) {
          if(AMPI_IN_PLACE != bufferSend) {
            datatype->copyIntoModifiedBuffer(bufferSend, 0, bufferSendMod, 0, count);
          } else {
            datatype->copyIntoModifiedBuffer(bufferRecv, 0, bufferRecvMod, 0, count);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(root == getCommRank(comm)) {
          if(AMPI_IN_PLACE != bufferSend) {
            h->bufferSendCount = datatype->computeActiveElements(count);
          } else {
            h->bufferSendCount = datatype->computeActiveElements(count);
          }
          h->bufferSendTotalSize = datatype->computeActiveElements(bufferSendElements);
          datatype->getADTool().createIndexTypeBuffer(h->bufferSendIndices, h->bufferSendTotalSize);
        }
        h->bufferRecvCount = datatype->computeActiveElements(count);
        h->bufferRecvTotalSize = datatype->computeActiveElements(bufferRecvElements);
        datatype->getADTool().createIndexTypeBuffer(h->bufferRecvIndices, h->bufferRecvTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          datatype->getADTool().createPrimalTypeBuffer(h->bufferRecvOldPrimals, h->bufferRecvTotalSize);
          datatype->getValues(bufferRecv, 0, h->bufferRecvOldPrimals, 0, count);
        }


        if(root == getCommRank(comm)) {
          if(AMPI_IN_PLACE != bufferSend) {
            datatype->getIndices(bufferSend, 0, h->bufferSendIndices, 0, count);
          } else {
            datatype->getIndices(bufferRecv, 0, h->bufferSendIndices, 0, count);
          }
        }

        datatype->createIndices(bufferRecv, 0, h->bufferRecvIndices, 0, count);

        // pack all the variables in the handle
        h->funcReverse = AMPI_Ibcast_wrap_b<DATATYPE>;
        h->funcForward = AMPI_Ibcast_wrap_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Ibcast_wrap_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->root = root;
        h->comm = comm;
      }

      if(!datatype->isModifiedBufferRequired()) {
        datatype->clearIndices(bufferRecv, 0, count);
      }

      rStatus = MPI_Ibcast_wrap(bufferSendMod, bufferRecvMod, count, datatype->getModifiedMpiType(), root, comm,
                                &request->request);

      AMPI_Ibcast_wrap_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Ibcast_wrap_AsyncHandle<DATATYPE>();
      asyncHandle->bufferSend = bufferSend;
      asyncHandle->bufferSendMod = bufferSendMod;
      asyncHandle->bufferRecv = bufferRecv;
      asyncHandle->bufferRecvMod = bufferRecvMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->root = root;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Ibcast_wrap_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Ibcast_wrap_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Ibcast_wrap_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Ibcast_wrap_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ibcast_wrap_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Ibcast_wrap_AsyncHandle<DATATYPE>*>(handle);
    typename DATATYPE::Type* bufferSend = asyncHandle->bufferSend;
    typename DATATYPE::ModifiedType* bufferSendMod = asyncHandle->bufferSendMod;
    typename DATATYPE::Type* bufferRecv = asyncHandle->bufferRecv;
    typename DATATYPE::ModifiedType* bufferRecvMod = asyncHandle->bufferRecvMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    int root = asyncHandle->root;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ibcast_wrap_AdjointHandle<DATATYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(bufferSend); // Unused generated to ignore warnings
    MEDI_UNUSED(bufferSendMod); // Unused generated to ignore warnings
    MEDI_UNUSED(bufferRecv); // Unused generated to ignore warnings
    MEDI_UNUSED(bufferRecvMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(root); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      datatype->getADTool().addToolAction(h);

      if(datatype->isModifiedBufferRequired()) {
        datatype->copyFromModifiedBuffer(bufferRecv, 0, bufferRecvMod, 0, count);
      }

      if(nullptr != h) {
        // handle the recv buffers
        datatype->registerValue(bufferRecv, 0, h->bufferRecvIndices, h->bufferRecvOldPrimals, 0, count);
      }

      datatype->getADTool().stopAssembly(h);

      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == bufferSend)) {
          datatype->deleteModifiedTypeBuffer(bufferSendMod);
        }
      }
      if(datatype->isModifiedBufferRequired() ) {
        datatype->deleteModifiedTypeBuffer(bufferRecvMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Igather_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Igather_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Igather_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    int sendcount;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Igather_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    if(h->root == getCommRank(h->comm)) {
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
      adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Igather_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Igather_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igather_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Igather_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Igather(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                   typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Igather(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), root,
                            comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcount;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        recvbufElements = recvcount * getCommSize(comm);

        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
        } else {
          recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
        }
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          recvtype->copyIntoModifiedBuffer(recvbuf, recvcount * getCommRank(comm), recvbufMod, recvcount * getCommRank(comm),
                                           recvcount);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(recvcount);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        if(root == getCommRank(comm)) {
          h->recvbufCount = recvtype->computeActiveElements(recvcount);
          h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);
        }


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          if(root == getCommRank(comm)) {
            recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
            if(root == getCommRank(comm)) {
              recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
            }
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          recvtype->getIndices(recvbuf, recvcount * getCommRank(comm), h->sendbufIndices, 0, recvcount);
        }

        if(root == getCommRank(comm)) {
          recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount * getCommSize(comm));
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Igather_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Igather_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Igather_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(root == getCommRank(comm)) {
        if(!recvtype->isModifiedBufferRequired()) {
          recvtype->clearIndices(recvbuf, 0, recvcount * getCommSize(comm));
        }
      }

      rStatus = MPI_Igather(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                            recvtype->getModifiedMpiType(), root, comm, &request->request);

      AMPI_Igather_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Igather_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sendcount = sendcount;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->recvcount = recvcount;
      asyncHandle->recvtype = recvtype;
      asyncHandle->root = root;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Igather_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Igather_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Igather_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Igather_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Igather_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = static_cast<AMPI_Igather_AsyncHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    int sendcount = asyncHandle->sendcount;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int recvcount = asyncHandle->recvcount;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    int root = asyncHandle->root;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igather_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcount); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcount); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(root); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired()) {
          recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount * getCommSize(comm));
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(root == getCommRank(comm)) {
          recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount * getCommSize(comm));
        }
      }

      recvtype->getADTool().stopAssembly(h);

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->deleteModifiedTypeBuffer(recvbufMod);
        }
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Igatherv_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int* recvbufCount;
    /* required for async */ int* recvbufCountVec;
    /* required for async */ int* recvbufDisplsVec;
    const  int* recvcounts;
    const  int* displs;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Igatherv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
      if(nullptr != recvbufCount) {
        delete [] recvbufCount;
        recvbufCount = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Igatherv_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    int sendcount;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    const int* displsMod;
    const  int* recvcounts;
    const  int* displs;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Igatherv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufPrimals, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->root, h->comm,
                                          &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    if(h->root == getCommRank(h->comm)) {
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
      adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
      delete [] h->recvbufCountVec;
      delete [] h->recvbufDisplsVec;
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Igatherv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->root, h->comm,
                                          &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
      delete [] h->recvbufCountVec;
      delete [] h->recvbufDisplsVec;
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->recvbufCountVec, h->recvbufDisplsVec, h->recvbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Igatherv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufAdjoints, h->recvbufCountVec, h->recvbufDisplsVec, h->recvcounts, h->displs, h->recvtype, h->root, h->comm,
                                          &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Igatherv_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
      delete [] h->recvbufCountVec;
      delete [] h->recvbufDisplsVec;
    }
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Igatherv_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Igatherv(MEDI_OPTIONAL_CONST typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype,
                    typename RECVTYPE::Type* recvbuf, const int* recvcounts, const int* displs, RECVTYPE* recvtype, int root,
                    AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Igatherv(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcounts, displs, recvtype->getMpiType(),
                             root, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* displsMod = displs;
      int displsTotalSize = 0;
      if(nullptr != displs) {
        displsTotalSize = computeDisplacementsTotalSize(recvcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          displsMod = createLinearDisplacements(recvcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = sendcount;
      } else {
        sendbufElements = recvcounts[getCommRank(comm)];
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        recvbufElements = displsTotalSize;

        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
        } else {
          recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
        }
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(sendtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->copyIntoModifiedBuffer(recvbuf, displs[rank], recvbufMod, displsMod[rank], recvcounts[rank]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
        } else {
          h->sendbufCount = recvtype->computeActiveElements(displs[getCommRank(comm)] + recvcounts[getCommRank(
                              comm)]) - recvtype->computeActiveElements(displs[getCommRank(comm)]);
        }
        h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        if(root == getCommRank(comm)) {
          createLinearIndexCounts(h->recvbufCount, recvcounts, displs, getCommSize(comm), recvtype);
          h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);
        }


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          if(root == getCommRank(comm)) {
            recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
            if(root == getCommRank(comm)) {
              for(int i = 0; i < getCommSize(comm); ++i) {
                recvtype->getValues(recvbuf, displs[i], h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
              }
            }
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount);
        } else {
          {
            const int rank = getCommRank(comm);
            recvtype->getIndices(recvbuf, displs[rank], h->sendbufIndices, 0, recvcounts[rank]);
          }
        }

        if(root == getCommRank(comm)) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->createIndices(recvbuf, displs[i], h->recvbufIndices, displsMod[i], recvcounts[i]);
          }
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Igatherv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Igatherv_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Igatherv_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcounts = recvcounts;
        h->displs = displs;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(root == getCommRank(comm)) {
        if(!recvtype->isModifiedBufferRequired()) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->clearIndices(recvbuf, displs[i], recvcounts[i]);
          }
        }
      }

      rStatus = MPI_Igatherv(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcounts, displsMod,
                             recvtype->getModifiedMpiType(), root, comm, &request->request);

      AMPI_Igatherv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Igatherv_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sendcount = sendcount;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->displsMod = displsMod;
      asyncHandle->recvcounts = recvcounts;
      asyncHandle->displs = displs;
      asyncHandle->recvtype = recvtype;
      asyncHandle->root = root;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Igatherv_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Igatherv_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Igatherv_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Igatherv_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Igatherv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Igatherv_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    MEDI_OPTIONAL_CONST  typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    int sendcount = asyncHandle->sendcount;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    const int* displsMod = asyncHandle->displsMod;
    const  int* recvcounts = asyncHandle->recvcounts;
    const  int* displs = asyncHandle->displs;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    int root = asyncHandle->root;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Igatherv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcount); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(displsMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcounts); // Unused generated to ignore warnings
    MEDI_UNUSED(displs); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(root); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired()) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->copyFromModifiedBuffer(recvbuf, displs[i], recvbufMod, displsMod[i], recvcounts[i]);
          }
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(root == getCommRank(comm)) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            recvtype->registerValue(recvbuf, displs[i], h->recvbufIndices, h->recvbufOldPrimals, displsMod[i], recvcounts[i]);
          }
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] displsMod;
      }

      if(sendtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        sendtype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(root == getCommRank(comm)) {
        if(recvtype->isModifiedBufferRequired() ) {
          recvtype->deleteModifiedTypeBuffer(recvbufMod);
        }
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Ireduce_global_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename DATATYPE::IndexType* sendbufIndices;
    typename DATATYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int recvbufTotalSize;
    typename DATATYPE::IndexType* recvbufIndices;
    typename DATATYPE::PrimalType* recvbufPrimals;
    typename DATATYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int count;
    DATATYPE* datatype;
    AMPI_Op op;
    int root;
    AMPI_Comm comm;

    ~AMPI_Ireduce_global_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename DATATYPE>
  struct AMPI_Ireduce_global_AsyncHandle : public AsyncHandle {
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* sendbuf;
    typename DATATYPE::ModifiedType* sendbufMod;
    typename DATATYPE::Type* recvbuf;
    typename DATATYPE::ModifiedType* recvbufMod;
    int count;
    DATATYPE* datatype;
    AMPI_Op op;
    int root;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename DATATYPE>
  void AMPI_Ireduce_global_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Ireduce_global_pri<DATATYPE>(h->sendbufPrimals, h->sendbufCountVec, h->recvbufPrimals, h->recvbufCountVec,
                                      h->count, h->datatype, h->op, h->root, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ireduce_global_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    if(h->root == getCommRank(h->comm)) {
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
      adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    }
  }

  template<typename DATATYPE>
  void AMPI_Ireduce_global_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Ireduce_global_fwd<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                      h->count, h->datatype, h->op, h->root, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ireduce_global_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename DATATYPE>
  void AMPI_Ireduce_global_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

      convOp.preAdjointOperation(h->recvbufAdjoints, h->recvbufPrimals, h->recvbufCount, adjointInterface->getVectorSize());
    }
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Ireduce_global_adj<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                      h->count, h->datatype, h->op, h->root, h->comm, &h->requestReverse);

  }

  template<typename DATATYPE>
  void AMPI_Ireduce_global_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>(handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    // the primals of the recive buffer are always given to the function. The operator should ignore them if not needed.
    // The wrapper functions make sure that for operators that need the primals an all* action is perfomed (e.g. Allreduce instead of Reduce)
    convOp.postAdjointOperation(h->sendbufAdjoints, h->sendbufPrimals, h->recvbufPrimals, h->sendbufTotalSize,
                                adjointInterface->getVectorSize());
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename DATATYPE>
  int AMPI_Ireduce_global_finish(HandleBase* handle);
  template<typename DATATYPE>
  int AMPI_Ireduce_global(MEDI_OPTIONAL_CONST typename DATATYPE::Type* sendbuf, typename DATATYPE::Type* recvbuf,
                          int count, DATATYPE* datatype, AMPI_Op op, int root, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;
    AMPI_Op convOp = datatype->getADTool().convertOperator(op);
    (void)convOp;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Ireduce(sendbuf, recvbuf, count, datatype->getMpiType(), convOp.primalFunction, root, comm,
                            &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = count;
      } else {
        sendbufElements = count;
      }

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(sendbuf));
      }
      typename DATATYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        recvbufElements = count;

        if(datatype->isModifiedBufferRequired() ) {
          datatype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
        } else {
          recvbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(recvbuf));
        }
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Ireduce_global_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          datatype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, count);
        } else {
          datatype->copyIntoModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = datatype->computeActiveElements(count);
        } else {
          h->sendbufCount = datatype->computeActiveElements(count);
        }
        h->sendbufTotalSize = datatype->computeActiveElements(sendbufElements);
        datatype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        if(root == getCommRank(comm)) {
          h->recvbufCount = datatype->computeActiveElements(count);
          h->recvbufTotalSize = datatype->computeActiveElements(recvbufElements);
          datatype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);
        }

        // extract the primal values for the operator if required
        if(convOp.requiresPrimal) {
          datatype->getADTool().createPrimalTypeBuffer(h->sendbufPrimals, h->sendbufTotalSize);
          if(AMPI_IN_PLACE != sendbuf) {
            datatype->getValues(sendbuf, 0, h->sendbufPrimals, 0, count);
          } else {
            datatype->getValues(recvbuf, 0, h->sendbufPrimals, 0, count);
          }
        }

        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          if(root == getCommRank(comm)) {
            datatype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
            if(root == getCommRank(comm)) {
              datatype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, count);
            }
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          datatype->getIndices(sendbuf, 0, h->sendbufIndices, 0, count);
        } else {
          datatype->getIndices(recvbuf, 0, h->sendbufIndices, 0, count);
        }

        if(root == getCommRank(comm)) {
          datatype->createIndices(recvbuf, 0, h->recvbufIndices, 0, count);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Ireduce_global_b<DATATYPE>;
        h->funcForward = AMPI_Ireduce_global_d_finish<DATATYPE>;
        h->funcPrimal = AMPI_Ireduce_global_p_finish<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->op = op;
        h->root = root;
        h->comm = comm;
      }

      if(root == getCommRank(comm)) {
        if(!datatype->isModifiedBufferRequired()) {
          datatype->clearIndices(recvbuf, 0, count);
        }
      }

      rStatus = MPI_Ireduce(sendbufMod, recvbufMod, count, datatype->getModifiedMpiType(), convOp.modifiedPrimalFunction,
                            root, comm, &request->request);

      AMPI_Ireduce_global_AsyncHandle<DATATYPE>* asyncHandle = new AMPI_Ireduce_global_AsyncHandle<DATATYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->count = count;
      asyncHandle->datatype = datatype;
      asyncHandle->op = op;
      asyncHandle->root = root;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Ireduce_global_finish<DATATYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Ireduce_global_b_finish<DATATYPE>,
                                           (ForwardFunction)AMPI_Ireduce_global_d<DATATYPE>, h);
        datatype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename DATATYPE>
  int AMPI_Ireduce_global_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Ireduce_global_AsyncHandle<DATATYPE>* asyncHandle = static_cast<AMPI_Ireduce_global_AsyncHandle<DATATYPE>*>
        (handle);
    MEDI_OPTIONAL_CONST  typename DATATYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename DATATYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    typename DATATYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename DATATYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int count = asyncHandle->count;
    DATATYPE* datatype = asyncHandle->datatype;
    AMPI_Op op = asyncHandle->op;
    int root = asyncHandle->root;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Ireduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Ireduce_global_AdjointHandle<DATATYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(count); // Unused generated to ignore warnings
    MEDI_UNUSED(datatype); // Unused generated to ignore warnings
    MEDI_UNUSED(op); // Unused generated to ignore warnings
    MEDI_UNUSED(root); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(datatype->getADTool().isActiveType()) {

      AMPI_Op convOp = datatype->getADTool().convertOperator(op);
      (void)convOp;
      datatype->getADTool().addToolAction(h);

      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired()) {
          datatype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(root == getCommRank(comm)) {
          datatype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, count);
        }
      }
      // extract the primal values for the operator if required
      if(nullptr != h && convOp.requiresPrimal) {
        if(root == getCommRank(comm)) {
          datatype->getADTool().createPrimalTypeBuffer(h->recvbufPrimals, h->recvbufTotalSize);
          if(root == getCommRank(comm)) {
            datatype->getValues(recvbuf, 0, h->recvbufPrimals, 0, count);
          }
        }
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired() ) {
          datatype->deleteModifiedTypeBuffer(recvbufMod);
        }
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iscatter_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Iscatter_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iscatter_AsyncHandle : public AsyncHandle {
    typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    int sendcount;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);

    }

    AMPI_Iscatter_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);

    }

    AMPI_Iscatter_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    }
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    }

    AMPI_Iscatter_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                          h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm, &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatter_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    }
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iscatter_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iscatter(typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::Type* recvbuf,
                    int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Iscatter(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), root,
                             comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        sendbufElements = sendcount * getCommSize(comm);

        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
        } else {
          sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
        }
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != recvbuf) {
        recvbufElements = recvcount;
      } else {
        recvbufElements = sendcount;
      }

      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired()) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount * getCommSize(comm));
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(root == getCommRank(comm)) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
          h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        }
        if(AMPI_IN_PLACE != recvbuf) {
          h->recvbufCount = recvtype->computeActiveElements(recvcount);
        } else {
          h->recvbufCount = sendtype->computeActiveElements(sendcount);
        }
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          if(AMPI_IN_PLACE != recvbuf) {
            recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount);
          } else {
            sendtype->getValues(sendbuf, sendcount * getCommRank(comm), h->recvbufOldPrimals, 0, sendcount);
          }
        }


        if(root == getCommRank(comm)) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount * getCommSize(comm));
        }

        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount);
        } else {
          sendtype->createIndices(sendbuf, sendcount * getCommRank(comm), h->recvbufIndices, 0, sendcount);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Iscatter_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Iscatter_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Iscatter_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->clearIndices(recvbuf, 0, recvcount);
        } else {
          sendtype->clearIndices(sendbuf, sendcount * getCommRank(comm), sendcount);
        }
      }

      rStatus = MPI_Iscatter(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                             recvtype->getModifiedMpiType(), root, comm, &request->request);

      AMPI_Iscatter_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Iscatter_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->sendcount = sendcount;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->recvcount = recvcount;
      asyncHandle->recvtype = recvtype;
      asyncHandle->root = root;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Iscatter_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Iscatter_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Iscatter_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iscatter_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Iscatter_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Iscatter_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    int sendcount = asyncHandle->sendcount;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int recvcount = asyncHandle->recvcount;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    int root = asyncHandle->root;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcount); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcount); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(root); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount);
        } else {
          sendtype->copyFromModifiedBuffer(sendbuf, sendcount * getCommRank(comm), sendbufMod, sendcount * getCommRank(comm),
                                           sendcount);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount);
        } else {
          sendtype->registerValue(sendbuf, sendcount * getCommRank(comm), h->recvbufIndices, h->recvbufOldPrimals, 0, sendcount);
        }
      }

      recvtype->getADTool().stopAssembly(h);

      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->deleteModifiedTypeBuffer(sendbufMod);
        }
      }
      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iscatterv_AdjointHandle : public AsyncAdjointHandle {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int* sendbufCount;
    /* required for async */ int* sendbufCountVec;
    /* required for async */ int* sendbufDisplsVec;
    const  int* sendcounts;
    const  int* displs;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Iscatterv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != sendbufCount) {
        delete [] sendbufCount;
        sendbufCount = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };

  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Iscatterv_AsyncHandle : public AsyncHandle {
    typename SENDTYPE::Type* sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod;
    const int* displsMod;
    const  int* sendcounts;
    const  int* displs;
    SENDTYPE* sendtype;
    typename RECVTYPE::Type* recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;
    AMPI_Request* request;
  };

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);

    }

    AMPI_Iscatterv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                           h->displs, h->sendtype, h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm,
                                           &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_p_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
      delete [] h->sendbufCountVec;
      delete [] h->sendbufDisplsVec;
    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);

    }

    AMPI_Iscatterv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                           h->displs, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm,
                                           &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_d_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
      delete [] h->sendbufCountVec;
      delete [] h->sendbufDisplsVec;
    }
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    }

    AMPI_Iscatterv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                           h->displs, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm,
                                           &h->requestReverse);

  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Iscatterv_b_finish(HandleBase* handle, AdjointInterface* adjointInterface) {

    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);
    MPI_Wait(&h->requestReverse.request, MPI_STATUS_IGNORE);

    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
      delete [] h->sendbufCountVec;
      delete [] h->sendbufDisplsVec;
    }
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iscatterv_finish(HandleBase* handle);
  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iscatterv(typename SENDTYPE::Type* sendbuf, const int* sendcounts, const int* displs, SENDTYPE* sendtype,
                     typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm, AMPI_Request* request) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Iscatterv(sendbuf, sendcounts, displs, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(),
                              root, comm, &request->request);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* displsMod = displs;
      int displsTotalSize = 0;
      if(nullptr != displs) {
        displsTotalSize = computeDisplacementsTotalSize(sendcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          displsMod = createLinearDisplacements(sendcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        sendbufElements = displsTotalSize;

        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
        } else {
          sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
        }
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != recvbuf) {
        recvbufElements = recvcount;
      } else {
        recvbufElements = sendcounts[getCommRank(comm)];
      }

      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired()) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->copyIntoModifiedBuffer(sendbuf, displs[i], sendbufMod, displsMod[i], sendcounts[i]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(root == getCommRank(comm)) {
          createLinearIndexCounts(h->sendbufCount, sendcounts, displs, getCommSize(comm), sendtype);
          h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        }
        if(AMPI_IN_PLACE != recvbuf) {
          h->recvbufCount = recvtype->computeActiveElements(recvcount);
        } else {
          h->recvbufCount = sendtype->computeActiveElements(displs[getCommRank(comm)] + sendcounts[getCommRank(
                              comm)]) - sendtype->computeActiveElements(displs[getCommRank(comm)]);
        }
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          if(AMPI_IN_PLACE != recvbuf) {
            recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount);
          } else {
            {
              const int rank = getCommRank(comm);
              sendtype->getValues(sendbuf, displs[rank], h->recvbufOldPrimals, 0, sendcounts[rank]);
            }
          }
        }


        if(root == getCommRank(comm)) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->getIndices(sendbuf, displs[i], h->sendbufIndices, displsMod[i], sendcounts[i]);
          }
        }

        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->createIndices(sendbuf, displs[rank], h->recvbufIndices, 0, sendcounts[rank]);
          }
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Iscatterv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Iscatterv_d_finish<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Iscatterv_p_finish<SENDTYPE, RECVTYPE>;
        h->sendcounts = sendcounts;
        h->displs = displs;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->clearIndices(recvbuf, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->clearIndices(sendbuf, displs[rank], sendcounts[rank]);
          }
        }
      }

      rStatus = MPI_Iscatterv(sendbufMod, sendcounts, displsMod, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                              recvtype->getModifiedMpiType(), root, comm, &request->request);

      AMPI_Iscatterv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle = new AMPI_Iscatterv_AsyncHandle<SENDTYPE, RECVTYPE>();
      asyncHandle->sendbuf = sendbuf;
      asyncHandle->sendbufMod = sendbufMod;
      asyncHandle->displsMod = displsMod;
      asyncHandle->sendcounts = sendcounts;
      asyncHandle->displs = displs;
      asyncHandle->sendtype = sendtype;
      asyncHandle->recvbuf = recvbuf;
      asyncHandle->recvbufMod = recvbufMod;
      asyncHandle->recvcount = recvcount;
      asyncHandle->recvtype = recvtype;
      asyncHandle->root = root;
      asyncHandle->comm = comm;
      asyncHandle->toolHandle = h;
      request->handle = asyncHandle;
      request->func = (ContinueFunction)AMPI_Iscatterv_finish<SENDTYPE, RECVTYPE>;

      // create adjoint wait
      if(nullptr != h) {
        WaitHandle* waitH = new WaitHandle((ReverseFunction)AMPI_Iscatterv_b_finish<SENDTYPE, RECVTYPE>,
                                           (ForwardFunction)AMPI_Iscatterv_d<SENDTYPE, RECVTYPE>, h);
        recvtype->getADTool().addToolAction(waitH);
      }
    }

    return rStatus;
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Iscatterv_finish(HandleBase* handle) {
    int rStatus = 0;

    AMPI_Iscatterv_AsyncHandle<SENDTYPE, RECVTYPE>* asyncHandle =
      static_cast<AMPI_Iscatterv_AsyncHandle<SENDTYPE, RECVTYPE>*>(handle);
    typename SENDTYPE::Type* sendbuf = asyncHandle->sendbuf;
    typename SENDTYPE::ModifiedType* sendbufMod = asyncHandle->sendbufMod;
    const int* displsMod = asyncHandle->displsMod;
    const  int* sendcounts = asyncHandle->sendcounts;
    const  int* displs = asyncHandle->displs;
    SENDTYPE* sendtype = asyncHandle->sendtype;
    typename RECVTYPE::Type* recvbuf = asyncHandle->recvbuf;
    typename RECVTYPE::ModifiedType* recvbufMod = asyncHandle->recvbufMod;
    int recvcount = asyncHandle->recvcount;
    RECVTYPE* recvtype = asyncHandle->recvtype;
    int root = asyncHandle->root;
    AMPI_Comm comm = asyncHandle->comm;
    AMPI_Request* request = asyncHandle->request;
    AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Iscatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (asyncHandle->toolHandle);
    MEDI_UNUSED(h); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(sendbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(displsMod); // Unused generated to ignore warnings
    MEDI_UNUSED(sendcounts); // Unused generated to ignore warnings
    MEDI_UNUSED(displs); // Unused generated to ignore warnings
    MEDI_UNUSED(sendtype); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbuf); // Unused generated to ignore warnings
    MEDI_UNUSED(recvbufMod); // Unused generated to ignore warnings
    MEDI_UNUSED(recvcount); // Unused generated to ignore warnings
    MEDI_UNUSED(recvtype); // Unused generated to ignore warnings
    MEDI_UNUSED(root); // Unused generated to ignore warnings
    MEDI_UNUSED(comm); // Unused generated to ignore warnings
    MEDI_UNUSED(request); // Unused generated to ignore warnings

    delete asyncHandle;

    if(recvtype->getADTool().isActiveType()) {

      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->copyFromModifiedBuffer(sendbuf, displs[rank], sendbufMod, displsMod[rank], sendcounts[rank]);
          }
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->registerValue(sendbuf, displs[rank], h->recvbufIndices, h->recvbufOldPrimals, 0, sendcounts[rank]);
          }
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] displsMod;
      }

      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->deleteModifiedTypeBuffer(sendbufMod);
        }
      }
      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  struct AMPI_Reduce_global_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename DATATYPE::IndexType* sendbufIndices;
    typename DATATYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int recvbufTotalSize;
    typename DATATYPE::IndexType* recvbufIndices;
    typename DATATYPE::PrimalType* recvbufPrimals;
    typename DATATYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int count;
    DATATYPE* datatype;
    AMPI_Op op;
    int root;
    AMPI_Comm comm;

    ~AMPI_Reduce_global_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        datatype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        datatype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename DATATYPE>
  void AMPI_Reduce_global_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Reduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Reduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);


    AMPI_Reduce_global_pri<DATATYPE>(h->sendbufPrimals, h->sendbufCountVec, h->recvbufPrimals, h->recvbufCountVec, h->count,
                                     h->datatype, h->op, h->root, h->comm);

    adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    if(h->root == getCommRank(h->comm)) {
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
      adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
    }
  }

  template<typename DATATYPE>
  void AMPI_Reduce_global_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Reduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Reduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);


    AMPI_Reduce_global_fwd<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                     h->count, h->datatype, h->op, h->root, h->comm);

    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename DATATYPE>
  void AMPI_Reduce_global_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Reduce_global_AdjointHandle<DATATYPE>* h = static_cast<AMPI_Reduce_global_AdjointHandle<DATATYPE>*>(handle);

    AMPI_Op convOp = h->datatype->getADTool().convertOperator(h->op);
    (void)convOp;
    h->recvbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
      adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

      convOp.preAdjointOperation(h->recvbufAdjoints, h->recvbufPrimals, h->recvbufCount, adjointInterface->getVectorSize());
    }
    if(h->datatype->getADTool().isOldPrimalsRequired()) {
      if(h->root == getCommRank(h->comm)) {
        adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
      }
    }
    h->sendbufAdjoints = nullptr;
    h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
    adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );

    AMPI_Reduce_global_adj<DATATYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->recvbufAdjoints, h->recvbufCountVec,
                                     h->count, h->datatype, h->op, h->root, h->comm);

    // the primals of the recive buffer are always given to the function. The operator should ignore them if not needed.
    // The wrapper functions make sure that for operators that need the primals an all* action is perfomed (e.g. Allreduce instead of Reduce)
    convOp.postAdjointOperation(h->sendbufAdjoints, h->sendbufPrimals, h->recvbufPrimals, h->sendbufTotalSize,
                                adjointInterface->getVectorSize());
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
    }
  }

  template<typename DATATYPE>
  int AMPI_Reduce_global(MEDI_OPTIONAL_CONST typename DATATYPE::Type* sendbuf, typename DATATYPE::Type* recvbuf,
                         int count, DATATYPE* datatype, AMPI_Op op, int root, AMPI_Comm comm) {
    int rStatus;
    AMPI_Op convOp = datatype->getADTool().convertOperator(op);
    (void)convOp;

    if(!datatype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Reduce(sendbuf, recvbuf, count, datatype->getMpiType(), convOp.primalFunction, root, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Reduce_global_AdjointHandle<DATATYPE>* h = nullptr;
      typename DATATYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != sendbuf) {
        sendbufElements = count;
      } else {
        sendbufElements = count;
      }

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
      } else {
        sendbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(sendbuf));
      }
      typename DATATYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        recvbufElements = count;

        if(datatype->isModifiedBufferRequired() ) {
          datatype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
        } else {
          recvbufMod = reinterpret_cast<typename DATATYPE::ModifiedType*>(const_cast<typename DATATYPE::Type*>(recvbuf));
        }
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(datatype->getADTool().isHandleRequired()) {
        h = new AMPI_Reduce_global_AdjointHandle<DATATYPE>();
      }
      datatype->getADTool().startAssembly(h);
      if(datatype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != sendbuf) {
          datatype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, count);
        } else {
          datatype->copyIntoModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(AMPI_IN_PLACE != sendbuf) {
          h->sendbufCount = datatype->computeActiveElements(count);
        } else {
          h->sendbufCount = datatype->computeActiveElements(count);
        }
        h->sendbufTotalSize = datatype->computeActiveElements(sendbufElements);
        datatype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        if(root == getCommRank(comm)) {
          h->recvbufCount = datatype->computeActiveElements(count);
          h->recvbufTotalSize = datatype->computeActiveElements(recvbufElements);
          datatype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);
        }

        // extract the primal values for the operator if required
        if(convOp.requiresPrimal) {
          datatype->getADTool().createPrimalTypeBuffer(h->sendbufPrimals, h->sendbufTotalSize);
          if(AMPI_IN_PLACE != sendbuf) {
            datatype->getValues(sendbuf, 0, h->sendbufPrimals, 0, count);
          } else {
            datatype->getValues(recvbuf, 0, h->sendbufPrimals, 0, count);
          }
        }

        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(datatype->getADTool().isOldPrimalsRequired()) {
          if(root == getCommRank(comm)) {
            datatype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
            if(root == getCommRank(comm)) {
              datatype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, count);
            }
          }
        }


        if(AMPI_IN_PLACE != sendbuf) {
          datatype->getIndices(sendbuf, 0, h->sendbufIndices, 0, count);
        } else {
          datatype->getIndices(recvbuf, 0, h->sendbufIndices, 0, count);
        }

        if(root == getCommRank(comm)) {
          datatype->createIndices(recvbuf, 0, h->recvbufIndices, 0, count);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Reduce_global_b<DATATYPE>;
        h->funcForward = AMPI_Reduce_global_d<DATATYPE>;
        h->funcPrimal = AMPI_Reduce_global_p<DATATYPE>;
        h->count = count;
        h->datatype = datatype;
        h->op = op;
        h->root = root;
        h->comm = comm;
      }

      if(root == getCommRank(comm)) {
        if(!datatype->isModifiedBufferRequired()) {
          datatype->clearIndices(recvbuf, 0, count);
        }
      }

      rStatus = MPI_Reduce(sendbufMod, recvbufMod, count, datatype->getModifiedMpiType(), convOp.modifiedPrimalFunction, root,
                           comm);
      datatype->getADTool().addToolAction(h);

      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired()) {
          datatype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, count);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(root == getCommRank(comm)) {
          datatype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, count);
        }
      }
      // extract the primal values for the operator if required
      if(nullptr != h && convOp.requiresPrimal) {
        if(root == getCommRank(comm)) {
          datatype->getADTool().createPrimalTypeBuffer(h->recvbufPrimals, h->recvbufTotalSize);
          if(root == getCommRank(comm)) {
            datatype->getValues(recvbuf, 0, h->recvbufPrimals, 0, count);
          }
        }
      }

      datatype->getADTool().stopAssembly(h);

      if(datatype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == sendbuf)) {
        datatype->deleteModifiedTypeBuffer(sendbufMod);
      }
      if(root == getCommRank(comm)) {
        if(datatype->isModifiedBufferRequired() ) {
          datatype->deleteModifiedTypeBuffer(recvbufMod);
        }
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Scatter_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int sendbufCount;
    int sendbufCountVec;
    int sendcount;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Scatter_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatter_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
      adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);

    }

    AMPI_Scatter_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatter_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);

    }

    AMPI_Scatter_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    }
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatter_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      h->sendbufCountVec = adjointInterface->getVectorSize() * h->sendbufCount;
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    }

    AMPI_Scatter_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendcount, h->sendtype,
                                         h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
    }
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Scatter(typename SENDTYPE::Type* sendbuf, int sendcount, SENDTYPE* sendtype, typename RECVTYPE::Type* recvbuf,
                   int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Scatter(sendbuf, sendcount, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(), root,
                            comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        sendbufElements = sendcount * getCommSize(comm);

        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
        } else {
          sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
        }
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != recvbuf) {
        recvbufElements = recvcount;
      } else {
        recvbufElements = sendcount;
      }

      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Scatter_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired()) {
          sendtype->copyIntoModifiedBuffer(sendbuf, 0, sendbufMod, 0, sendcount * getCommSize(comm));
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(root == getCommRank(comm)) {
          h->sendbufCount = sendtype->computeActiveElements(sendcount);
          h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        }
        if(AMPI_IN_PLACE != recvbuf) {
          h->recvbufCount = recvtype->computeActiveElements(recvcount);
        } else {
          h->recvbufCount = sendtype->computeActiveElements(sendcount);
        }
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          if(AMPI_IN_PLACE != recvbuf) {
            recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount);
          } else {
            sendtype->getValues(sendbuf, sendcount * getCommRank(comm), h->recvbufOldPrimals, 0, sendcount);
          }
        }


        if(root == getCommRank(comm)) {
          sendtype->getIndices(sendbuf, 0, h->sendbufIndices, 0, sendcount * getCommSize(comm));
        }

        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount);
        } else {
          sendtype->createIndices(sendbuf, sendcount * getCommRank(comm), h->recvbufIndices, 0, sendcount);
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Scatter_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Scatter_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Scatter_p<SENDTYPE, RECVTYPE>;
        h->sendcount = sendcount;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->clearIndices(recvbuf, 0, recvcount);
        } else {
          sendtype->clearIndices(sendbuf, sendcount * getCommRank(comm), sendcount);
        }
      }

      rStatus = MPI_Scatter(sendbufMod, sendcount, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                            recvtype->getModifiedMpiType(), root, comm);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount);
        } else {
          sendtype->copyFromModifiedBuffer(sendbuf, sendcount * getCommRank(comm), sendbufMod, sendcount * getCommRank(comm),
                                           sendcount);
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount);
        } else {
          sendtype->registerValue(sendbuf, sendcount * getCommRank(comm), h->recvbufIndices, h->recvbufOldPrimals, 0, sendcount);
        }
      }

      recvtype->getADTool().stopAssembly(h);

      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->deleteModifiedTypeBuffer(sendbufMod);
        }
      }
      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename SENDTYPE, typename RECVTYPE>
  struct AMPI_Scatterv_AdjointHandle : public HandleBase {
    int sendbufTotalSize;
    typename SENDTYPE::IndexType* sendbufIndices;
    typename SENDTYPE::PrimalType* sendbufPrimals;
    /* required for async */ void* sendbufAdjoints;
    int* sendbufCount;
    /* required for async */ int* sendbufCountVec;
    /* required for async */ int* sendbufDisplsVec;
    const  int* sendcounts;
    const  int* displs;
    SENDTYPE* sendtype;
    int recvbufTotalSize;
    typename RECVTYPE::IndexType* recvbufIndices;
    typename RECVTYPE::PrimalType* recvbufPrimals;
    typename RECVTYPE::PrimalType* recvbufOldPrimals;
    /* required for async */ void* recvbufAdjoints;
    int recvbufCount;
    int recvbufCountVec;
    int recvcount;
    RECVTYPE* recvtype;
    int root;
    AMPI_Comm comm;

    ~AMPI_Scatterv_AdjointHandle () {
      if(nullptr != sendbufIndices) {
        sendtype->getADTool().deleteIndexTypeBuffer(sendbufIndices);
        sendbufIndices = nullptr;
      }
      if(nullptr != sendbufPrimals) {
        sendtype->getADTool().deletePrimalTypeBuffer(sendbufPrimals);
        sendbufPrimals = nullptr;
      }
      if(nullptr != sendbufCount) {
        delete [] sendbufCount;
        sendbufCount = nullptr;
      }
      if(nullptr != recvbufIndices) {
        recvtype->getADTool().deleteIndexTypeBuffer(recvbufIndices);
        recvbufIndices = nullptr;
      }
      if(nullptr != recvbufPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufPrimals);
        recvbufPrimals = nullptr;
      }
      if(nullptr != recvbufOldPrimals) {
        recvtype->getADTool().deletePrimalTypeBuffer(recvbufOldPrimals);
        recvbufOldPrimals = nullptr;
      }
    }
  };


  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatterv_p(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createPrimalTypeBuffer((void*&)h->recvbufPrimals, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createPrimalTypeBuffer((void*&)h->sendbufPrimals, h->sendbufTotalSize );
      // Primal buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getPrimals(h->sendbufIndices, h->sendbufPrimals, h->sendbufTotalSize);

    }

    AMPI_Scatterv_pri<SENDTYPE, RECVTYPE>(h->sendbufPrimals, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                          h->displs, h->sendtype, h->recvbufPrimals, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deletePrimalTypeBuffer((void*&)h->sendbufPrimals);
      delete [] h->sendbufCountVec;
      delete [] h->sendbufDisplsVec;
    }
    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->getPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    // Primal buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->setPrimals(h->recvbufIndices, h->recvbufPrimals, h->recvbufTotalSize);
    adjointInterface->deletePrimalTypeBuffer((void*&)h->recvbufPrimals);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatterv_d(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->getAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);

    }

    AMPI_Scatterv_fwd<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                          h->displs, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
      delete [] h->sendbufCountVec;
      delete [] h->sendbufDisplsVec;
    }
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->updateAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  void AMPI_Scatterv_b(HandleBase* handle, AdjointInterface* adjointInterface) {
    AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = static_cast<AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>*>
        (handle);

    h->recvbufAdjoints = nullptr;
    h->recvbufCountVec = adjointInterface->getVectorSize() * h->recvbufCount;
    adjointInterface->createAdjointTypeBuffer(h->recvbufAdjoints, h->recvbufTotalSize );
    // Adjoint buffers are always linear in space so we can accesses them in one sweep
    adjointInterface->getAdjoints(h->recvbufIndices, h->recvbufAdjoints, h->recvbufTotalSize);

    if(h->recvtype->getADTool().isOldPrimalsRequired()) {
      adjointInterface->setPrimals(h->recvbufIndices, h->recvbufOldPrimals, h->recvbufTotalSize);
    }
    h->sendbufAdjoints = nullptr;
    if(h->root == getCommRank(h->comm)) {
      createLinearDisplacementsAndCount(h->sendbufCountVec, h->sendbufDisplsVec, h->sendbufCount, getCommSize(h->comm),
                                        adjointInterface->getVectorSize());
      adjointInterface->createAdjointTypeBuffer(h->sendbufAdjoints, h->sendbufTotalSize );
    }

    AMPI_Scatterv_adj<SENDTYPE, RECVTYPE>(h->sendbufAdjoints, h->sendbufCountVec, h->sendbufDisplsVec, h->sendcounts,
                                          h->displs, h->sendtype, h->recvbufAdjoints, h->recvbufCountVec, h->recvcount, h->recvtype, h->root, h->comm);

    if(h->root == getCommRank(h->comm)) {
      // Adjoint buffers are always linear in space so we can accesses them in one sweep
      adjointInterface->updateAdjoints(h->sendbufIndices, h->sendbufAdjoints, h->sendbufTotalSize);
      adjointInterface->deleteAdjointTypeBuffer(h->sendbufAdjoints);
      delete [] h->sendbufCountVec;
      delete [] h->sendbufDisplsVec;
    }
    adjointInterface->deleteAdjointTypeBuffer(h->recvbufAdjoints);
  }

  template<typename SENDTYPE, typename RECVTYPE>
  int AMPI_Scatterv(typename SENDTYPE::Type* sendbuf, const int* sendcounts, const int* displs, SENDTYPE* sendtype,
                    typename RECVTYPE::Type* recvbuf, int recvcount, RECVTYPE* recvtype, int root, AMPI_Comm comm) {
    int rStatus;

    if(!recvtype->getADTool().isActiveType()) {
      // call the regular function if the type is not active
      rStatus = MPI_Scatterv(sendbuf, sendcounts, displs, sendtype->getMpiType(), recvbuf, recvcount, recvtype->getMpiType(),
                             root, comm);
    } else {

      // the type is an AD type so handle the buffers
      AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>* h = nullptr;
      MEDI_OPTIONAL_CONST int* displsMod = displs;
      int displsTotalSize = 0;
      if(nullptr != displs) {
        displsTotalSize = computeDisplacementsTotalSize(sendcounts, getCommSize(comm));
        if(recvtype->isModifiedBufferRequired()) {
          displsMod = createLinearDisplacements(sendcounts, getCommSize(comm));
        }
      }
      typename SENDTYPE::ModifiedType* sendbufMod = nullptr;
      int sendbufElements = 0;

      if(root == getCommRank(comm)) {
        // compute the total size of the buffer
        sendbufElements = displsTotalSize;

        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->createModifiedTypeBuffer(sendbufMod, sendbufElements);
        } else {
          sendbufMod = reinterpret_cast<typename SENDTYPE::ModifiedType*>(const_cast<typename SENDTYPE::Type*>(sendbuf));
        }
      }
      typename RECVTYPE::ModifiedType* recvbufMod = nullptr;
      int recvbufElements = 0;

      // compute the total size of the buffer
      if(AMPI_IN_PLACE != recvbuf) {
        recvbufElements = recvcount;
      } else {
        recvbufElements = sendcounts[getCommRank(comm)];
      }

      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->createModifiedTypeBuffer(recvbufMod, recvbufElements);
      } else {
        recvbufMod = reinterpret_cast<typename RECVTYPE::ModifiedType*>(const_cast<typename RECVTYPE::Type*>(recvbuf));
      }

      // the handle is created if a reverse action should be recorded, h != nullptr => tape is active
      if(recvtype->getADTool().isHandleRequired()) {
        h = new AMPI_Scatterv_AdjointHandle<SENDTYPE, RECVTYPE>();
      }
      recvtype->getADTool().startAssembly(h);
      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired()) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->copyIntoModifiedBuffer(sendbuf, displs[i], sendbufMod, displsMod[i], sendcounts[i]);
          }
        }
      }

      if(nullptr != h) {
        // gather the information for the reverse sweep

        // create the index buffers
        if(root == getCommRank(comm)) {
          createLinearIndexCounts(h->sendbufCount, sendcounts, displs, getCommSize(comm), sendtype);
          h->sendbufTotalSize = sendtype->computeActiveElements(sendbufElements);
          recvtype->getADTool().createIndexTypeBuffer(h->sendbufIndices, h->sendbufTotalSize);
        }
        if(AMPI_IN_PLACE != recvbuf) {
          h->recvbufCount = recvtype->computeActiveElements(recvcount);
        } else {
          h->recvbufCount = sendtype->computeActiveElements(displs[getCommRank(comm)] + sendcounts[getCommRank(
                              comm)]) - sendtype->computeActiveElements(displs[getCommRank(comm)]);
        }
        h->recvbufTotalSize = recvtype->computeActiveElements(recvbufElements);
        recvtype->getADTool().createIndexTypeBuffer(h->recvbufIndices, h->recvbufTotalSize);


        // extract the old primal values from the recv buffer if the AD tool
        // needs the primal values reset
        if(recvtype->getADTool().isOldPrimalsRequired()) {
          recvtype->getADTool().createPrimalTypeBuffer(h->recvbufOldPrimals, h->recvbufTotalSize);
          if(AMPI_IN_PLACE != recvbuf) {
            recvtype->getValues(recvbuf, 0, h->recvbufOldPrimals, 0, recvcount);
          } else {
            {
              const int rank = getCommRank(comm);
              sendtype->getValues(sendbuf, displs[rank], h->recvbufOldPrimals, 0, sendcounts[rank]);
            }
          }
        }


        if(root == getCommRank(comm)) {
          for(int i = 0; i < getCommSize(comm); ++i) {
            sendtype->getIndices(sendbuf, displs[i], h->sendbufIndices, displsMod[i], sendcounts[i]);
          }
        }

        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->createIndices(recvbuf, 0, h->recvbufIndices, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->createIndices(sendbuf, displs[rank], h->recvbufIndices, 0, sendcounts[rank]);
          }
        }

        // pack all the variables in the handle
        h->funcReverse = AMPI_Scatterv_b<SENDTYPE, RECVTYPE>;
        h->funcForward = AMPI_Scatterv_d<SENDTYPE, RECVTYPE>;
        h->funcPrimal = AMPI_Scatterv_p<SENDTYPE, RECVTYPE>;
        h->sendcounts = sendcounts;
        h->displs = displs;
        h->sendtype = sendtype;
        h->recvcount = recvcount;
        h->recvtype = recvtype;
        h->root = root;
        h->comm = comm;
      }

      if(!recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->clearIndices(recvbuf, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->clearIndices(sendbuf, displs[rank], sendcounts[rank]);
          }
        }
      }

      rStatus = MPI_Scatterv(sendbufMod, sendcounts, displsMod, sendtype->getModifiedMpiType(), recvbufMod, recvcount,
                             recvtype->getModifiedMpiType(), root, comm);
      recvtype->getADTool().addToolAction(h);

      if(recvtype->isModifiedBufferRequired()) {
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->copyFromModifiedBuffer(recvbuf, 0, recvbufMod, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->copyFromModifiedBuffer(sendbuf, displs[rank], sendbufMod, displsMod[rank], sendcounts[rank]);
          }
        }
      }

      if(nullptr != h) {
        // handle the recv buffers
        if(AMPI_IN_PLACE != recvbuf) {
          recvtype->registerValue(recvbuf, 0, h->recvbufIndices, h->recvbufOldPrimals, 0, recvcount);
        } else {
          {
            const int rank = getCommRank(comm);
            sendtype->registerValue(sendbuf, displs[rank], h->recvbufIndices, h->recvbufOldPrimals, 0, sendcounts[rank]);
          }
        }
      }

      recvtype->getADTool().stopAssembly(h);
      if(recvtype->isModifiedBufferRequired()) {
        delete [] displsMod;
      }

      if(root == getCommRank(comm)) {
        if(sendtype->isModifiedBufferRequired() ) {
          sendtype->deleteModifiedTypeBuffer(sendbufMod);
        }
      }
      if(recvtype->isModifiedBufferRequired()  && !(AMPI_IN_PLACE == recvbuf)) {
        recvtype->deleteModifiedTypeBuffer(recvbufMod);
      }

      // handle is deleted by the AD tool
    }

    return rStatus;
  }

#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Buffer_attach(void* buffer, int size) {
    return MPI_Buffer_attach(buffer, size);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Buffer_detach(void* buffer_addr, int* size) {
    return MPI_Buffer_detach(buffer_addr, size);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cancel(AMPI_Request* request) {
    return MPI_Cancel(&request->request);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Get_count(MEDI_OPTIONAL_CONST AMPI_Status* status, DATATYPE* datatype, int* count) {
    return MPI_Get_count(status, datatype->getModifiedMpiType(), count);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Iprobe(int source, int tag, AMPI_Comm comm, int* flag, AMPI_Status* status) {
    return MPI_Iprobe(source, tag, comm, flag, status);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Probe(int source, int tag, AMPI_Comm comm, AMPI_Status* status) {
    return MPI_Probe(source, tag, comm, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Request_get_status(AMPI_Request request, int* flag, AMPI_Status* status) {
    return MPI_Request_get_status(request.request, flag, status);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Test_cancelled(MEDI_OPTIONAL_CONST AMPI_Status* status, int* flag) {
    return MPI_Test_cancelled(status, flag);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  inline MPI_Aint AMPI_Aint_add(AMPI_Aint base, AMPI_Aint disp) {
    return MPI_Aint_add(base, disp);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  inline MPI_Aint AMPI_Aint_diff(AMPI_Aint addr1, AMPI_Aint addr2) {
    return MPI_Aint_diff(addr1, addr2);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Get_address(MEDI_OPTIONAL_CONST void* location, AMPI_Aint* address) {
    return MPI_Get_address(location, address);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Get_elements(const AMPI_Status* status, DATATYPE* datatype, int* count) {
    return MPI_Get_elements(status, datatype->getModifiedMpiType(), count);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Get_elements_x(const AMPI_Status* status, DATATYPE* datatype, AMPI_Count* count) {
    return MPI_Get_elements_x(status, datatype->getModifiedMpiType(), count);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_contents(DATATYPE* datatype, int max_integers, int max_addresses, int max_datatypes,
                                    int* array_of_integers, AMPI_Aint* array_of_addresses, AMPI_Datatype* array_of_datatypes) {
    return MPI_Type_get_contents(datatype->getModifiedMpiType(), max_integers, max_addresses, max_datatypes,
                                 array_of_integers, array_of_addresses, array_of_datatypes);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_envelope(DATATYPE* datatype, int* num_integers, int* num_addresses, int* num_datatypes,
                                    int* combiner) {
    return MPI_Type_get_envelope(datatype->getModifiedMpiType(), num_integers, num_addresses, num_datatypes, combiner);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_extent(DATATYPE* datatype, AMPI_Aint* lb, AMPI_Aint* extent) {
    return MPI_Type_get_extent(datatype->getModifiedMpiType(), lb, extent);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_extent_x(DATATYPE* datatype, AMPI_Count* lb, AMPI_Count* extent) {
    return MPI_Type_get_extent_x(datatype->getModifiedMpiType(), lb, extent);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_true_extent(DATATYPE* datatype, AMPI_Aint* true_lb, AMPI_Aint* true_extent) {
    return MPI_Type_get_true_extent(datatype->getModifiedMpiType(), true_lb, true_extent);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_true_extent_x(DATATYPE* datatype, AMPI_Count* true_lb, AMPI_Count* true_extent) {
    return MPI_Type_get_true_extent_x(datatype->getModifiedMpiType(), true_lb, true_extent);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename OLDTYPE, typename NEWTYPE>
  inline int AMPI_Type_hindexed(int count, const int* array_of_blocklengths, const AMPI_Aint* array_of_displacements,
                                OLDTYPE* oldtype, NEWTYPE* newtype) {
    return MPI_Type_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype->getModifiedMpiType(),
                             newtype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename OLDTYPE, typename NEWTYPE>
  inline int AMPI_Type_hvector(int count, int blocklength, AMPI_Aint stride, OLDTYPE* oldtype, NEWTYPE* newtype) {
    return MPI_Type_hvector(count, blocklength, stride, oldtype->getModifiedMpiType(), newtype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_size(DATATYPE* datatype, int* size) {
    return MPI_Type_size(datatype->getModifiedMpiType(), size);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_size_x(DATATYPE* datatype, AMPI_Count* size) {
    return MPI_Type_size_x(datatype->getModifiedMpiType(), size);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  template<typename ARRAY_OF_TYPES, typename NEWTYPE>
  inline int AMPI_Type_struct(int count, const int* array_of_blocklengths, const AMPI_Aint* array_of_displacements,
                              const ARRAY_OF_TYPES* array_of_types, NEWTYPE* newtype) {
    return MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types->getModifiedMpiType(),
                           newtype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Barrier(AMPI_Comm comm) {
    return MPI_Barrier(comm);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Op_commutative(AMPI_Op op, int* commute) {
    AMPI_Op convOp = op;
    (void)convOp;
    return MPI_Op_commutative(convOp.modifiedPrimalFunction, commute);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_compare(AMPI_Comm comm1, AMPI_Comm comm2, int* result) {
    return MPI_Comm_compare(comm1, comm2, result);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_create(AMPI_Comm comm, AMPI_Group group, AMPI_Comm* newcomm) {
    return MPI_Comm_create(comm, group, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_create_group(AMPI_Comm comm, AMPI_Group group, int tag, AMPI_Comm* newcomm) {
    return MPI_Comm_create_group(comm, group, tag, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_create_keyval(AMPI_Comm_copy_attr_function* comm_copy_attr_fn,
                                     AMPI_Comm_delete_attr_function* comm_delete_attr_fn, int* comm_keyval, void* extra_state) {
    return MPI_Comm_create_keyval(comm_copy_attr_fn, comm_delete_attr_fn, comm_keyval, extra_state);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_delete_attr(AMPI_Comm comm, int comm_keyval) {
    return MPI_Comm_delete_attr(comm, comm_keyval);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_dup(AMPI_Comm comm, AMPI_Comm* newcomm) {
    return MPI_Comm_dup(comm, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_dup_with_info(AMPI_Comm comm, AMPI_Info info, AMPI_Comm* newcomm) {
    return MPI_Comm_dup_with_info(comm, info, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_free(AMPI_Comm* comm) {
    return MPI_Comm_free(comm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_free_keyval(int* comm_keyval) {
    return MPI_Comm_free_keyval(comm_keyval);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_get_attr(AMPI_Comm comm, int comm_keyval, void* attribute_val, int* flag) {
    return MPI_Comm_get_attr(comm, comm_keyval, attribute_val, flag);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_get_info(AMPI_Comm comm, AMPI_Info* info_used) {
    return MPI_Comm_get_info(comm, info_used);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_get_name(AMPI_Comm comm, char* comm_name, int* resultlen) {
    return MPI_Comm_get_name(comm, comm_name, resultlen);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_group(AMPI_Comm comm, AMPI_Group* group) {
    return MPI_Comm_group(comm, group);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_idup(AMPI_Comm comm, AMPI_Comm* newcomm, AMPI_Request* request) {
    return MPI_Comm_idup(comm, newcomm, &request->request);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_rank(AMPI_Comm comm, int* rank) {
    return MPI_Comm_rank(comm, rank);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_remote_group(AMPI_Comm comm, AMPI_Group* group) {
    return MPI_Comm_remote_group(comm, group);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_remote_size(AMPI_Comm comm, int* size) {
    return MPI_Comm_remote_size(comm, size);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_set_attr(AMPI_Comm comm, int comm_keyval, void* attribute_val) {
    return MPI_Comm_set_attr(comm, comm_keyval, attribute_val);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_set_info(AMPI_Comm comm, AMPI_Info info) {
    return MPI_Comm_set_info(comm, info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_set_name(AMPI_Comm comm, MEDI_OPTIONAL_CONST char* comm_name) {
    return MPI_Comm_set_name(comm, comm_name);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_size(AMPI_Comm comm, int* size) {
    return MPI_Comm_size(comm, size);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_split(AMPI_Comm comm, int color, int key, AMPI_Comm* newcomm) {
    return MPI_Comm_split(comm, color, key, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_split_type(AMPI_Comm comm, int split_type, int key, AMPI_Info info, AMPI_Comm* newcomm) {
    return MPI_Comm_split_type(comm, split_type, key, info, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_test_inter(AMPI_Comm comm, int* flag) {
    return MPI_Comm_test_inter(comm, flag);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_compare(AMPI_Group group1, AMPI_Group group2, int* result) {
    return MPI_Group_compare(group1, group2, result);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_difference(AMPI_Group group1, AMPI_Group group2, AMPI_Group* newgroup) {
    return MPI_Group_difference(group1, group2, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_excl(AMPI_Group group, int n, MEDI_OPTIONAL_CONST int* ranks, AMPI_Group* newgroup) {
    return MPI_Group_excl(group, n, ranks, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_free(AMPI_Group* group) {
    return MPI_Group_free(group);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_incl(AMPI_Group group, int n, MEDI_OPTIONAL_CONST int* ranks, AMPI_Group* newgroup) {
    return MPI_Group_incl(group, n, ranks, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_intersection(AMPI_Group group1, AMPI_Group group2, AMPI_Group* newgroup) {
    return MPI_Group_intersection(group1, group2, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_range_excl(AMPI_Group group, int n, ::medi::common::Range* ranges, AMPI_Group* newgroup) {
    return MPI_Group_range_excl(group, n, ranges, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_range_incl(AMPI_Group group, int n, ::medi::common::Range* ranges, AMPI_Group* newgroup) {
    return MPI_Group_range_incl(group, n, ranges, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_rank(AMPI_Group group, int* rank) {
    return MPI_Group_rank(group, rank);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_size(AMPI_Group group, int* size) {
    return MPI_Group_size(group, size);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_translate_ranks(AMPI_Group group1, int n, MEDI_OPTIONAL_CONST int* ranks1, AMPI_Group group2,
                                        int* ranks2) {
    return MPI_Group_translate_ranks(group1, n, ranks1, group2, ranks2);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Group_union(AMPI_Group group1, AMPI_Group group2, AMPI_Group* newgroup) {
    return MPI_Group_union(group1, group2, newgroup);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Intercomm_create(AMPI_Comm local_comm, int local_leader, AMPI_Comm peer_comm, int remote_leader,
                                   int tag, AMPI_Comm* newintercomm) {
    return MPI_Intercomm_create(local_comm, local_leader, peer_comm, remote_leader, tag, newintercomm);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Intercomm_merge(AMPI_Comm intercomm, int high, AMPI_Comm* newintracomm) {
    return MPI_Intercomm_merge(intercomm, high, newintracomm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Type_create_keyval(AMPI_Type_copy_attr_function* type_copy_attr_fn,
                                     AMPI_Type_delete_attr_function* type_delete_attr_fn, int* type_keyval, void* extra_state) {
    return MPI_Type_create_keyval(type_copy_attr_fn, type_delete_attr_fn, type_keyval, extra_state);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_delete_attr(DATATYPE* datatype, int type_keyval) {
    return MPI_Type_delete_attr(datatype->getModifiedMpiType(), type_keyval);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Type_free_keyval(int* type_keyval) {
    return MPI_Type_free_keyval(type_keyval);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_attr(DATATYPE* datatype, int type_keyval, void* attribute_val, int* flag) {
    return MPI_Type_get_attr(datatype->getModifiedMpiType(), type_keyval, attribute_val, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_get_name(DATATYPE* datatype, char* type_name, int* resultlen) {
    return MPI_Type_get_name(datatype->getModifiedMpiType(), type_name, resultlen);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_set_attr(DATATYPE* datatype, int type_keyval, void* attribute_val) {
    return MPI_Type_set_attr(datatype->getModifiedMpiType(), type_keyval, attribute_val);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Type_set_name(DATATYPE* datatype, const char* type_name) {
    return MPI_Type_set_name(datatype->getModifiedMpiType(), type_name);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_create_keyval(AMPI_Win_copy_attr_function* win_copy_attr_fn,
                                    AMPI_Win_delete_attr_function* win_delete_attr_fn, int* win_keyval, void* extra_state) {
    return MPI_Win_create_keyval(win_copy_attr_fn, win_delete_attr_fn, win_keyval, extra_state);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_delete_attr(AMPI_Win win, int win_keyval) {
    return MPI_Win_delete_attr(win, win_keyval);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_free_keyval(int* win_keyval) {
    return MPI_Win_free_keyval(win_keyval);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_get_attr(AMPI_Win win, int win_keyval, void* attribute_val, int* flag) {
    return MPI_Win_get_attr(win, win_keyval, attribute_val, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_get_name(AMPI_Win win, char* win_name, int* resultlen) {
    return MPI_Win_get_name(win, win_name, resultlen);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_set_attr(AMPI_Win win, int win_keyval, void* attribute_val) {
    return MPI_Win_set_attr(win, win_keyval, attribute_val);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_set_name(AMPI_Win win, MEDI_OPTIONAL_CONST char* win_name) {
    return MPI_Win_set_name(win, win_name);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_coords(AMPI_Comm comm, int rank, int maxdims, int* coords) {
    return MPI_Cart_coords(comm, rank, maxdims, coords);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_create(AMPI_Comm comm_old, int ndims, MEDI_OPTIONAL_CONST int* dims,
                              MEDI_OPTIONAL_CONST int* periods, int reorder, AMPI_Comm* comm_cart) {
    return MPI_Cart_create(comm_old, ndims, dims, periods, reorder, comm_cart);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_get(AMPI_Comm comm, int maxdims, int* dims, int* periods, int* coords) {
    return MPI_Cart_get(comm, maxdims, dims, periods, coords);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_map(AMPI_Comm comm, int ndims, MEDI_OPTIONAL_CONST int* dims, MEDI_OPTIONAL_CONST int* periods,
                           int* newrank) {
    return MPI_Cart_map(comm, ndims, dims, periods, newrank);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_rank(AMPI_Comm comm, MEDI_OPTIONAL_CONST int* coords, int* rank) {
    return MPI_Cart_rank(comm, coords, rank);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_shift(AMPI_Comm comm, int direction, int disp, int* rank_source, int* rank_dest) {
    return MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cart_sub(AMPI_Comm comm, MEDI_OPTIONAL_CONST int* remain_dims, AMPI_Comm* newcomm) {
    return MPI_Cart_sub(comm, remain_dims, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Cartdim_get(AMPI_Comm comm, int* ndims) {
    return MPI_Cartdim_get(comm, ndims);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Dims_create(int nnodes, int ndims, int* dims) {
    return MPI_Dims_create(nnodes, ndims, dims);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Dist_graph_create(AMPI_Comm comm_old, int n, MEDI_OPTIONAL_CONST int* sources,
                                    MEDI_OPTIONAL_CONST int* degrees, MEDI_OPTIONAL_CONST int* destinations, MEDI_OPTIONAL_CONST int* weights,
                                    AMPI_Info info, int reorder, AMPI_Comm* comm_dist_graph) {
    return MPI_Dist_graph_create(comm_old, n, sources, degrees, destinations, weights, info, reorder, comm_dist_graph);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Dist_graph_create_adjacent(AMPI_Comm comm_old, int indegree, MEDI_OPTIONAL_CONST int* sources,
      MEDI_OPTIONAL_CONST int* sourceweights, int outdegree, MEDI_OPTIONAL_CONST int* destinations,
      MEDI_OPTIONAL_CONST int* destweights, AMPI_Info info, int reorder, AMPI_Comm* comm_dist_graph) {
    return MPI_Dist_graph_create_adjacent(comm_old, indegree, sources, sourceweights, outdegree, destinations, destweights,
                                          info, reorder, comm_dist_graph);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Dist_graph_neighbors(AMPI_Comm comm, int maxindegree, int* sources, int* sourceweights,
                                       int maxoutdegree, int* destinations, int* destweights) {
    return MPI_Dist_graph_neighbors(comm, maxindegree, sources, sourceweights, maxoutdegree, destinations, destweights);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Dist_graph_neighbors_count(AMPI_Comm comm, int* indegree, int* outdegree, int* weighted) {
    return MPI_Dist_graph_neighbors_count(comm, indegree, outdegree, weighted);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Graph_create(AMPI_Comm comm_old, int nnodes, MEDI_OPTIONAL_CONST int* index,
                               MEDI_OPTIONAL_CONST int* edges, int reorder, AMPI_Comm* comm_graph) {
    return MPI_Graph_create(comm_old, nnodes, index, edges, reorder, comm_graph);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Graph_get(AMPI_Comm comm, int maxindex, int maxedges, int* index, int* edges) {
    return MPI_Graph_get(comm, maxindex, maxedges, index, edges);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Graph_map(AMPI_Comm comm, int nnodes, MEDI_OPTIONAL_CONST int* index, MEDI_OPTIONAL_CONST int* edges,
                            int* newrank) {
    return MPI_Graph_map(comm, nnodes, index, edges, newrank);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Graph_neighbors(AMPI_Comm comm, int rank, int maxneighbors, int* neighbors) {
    return MPI_Graph_neighbors(comm, rank, maxneighbors, neighbors);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Graph_neighbors_count(AMPI_Comm comm, int rank, int* nneighbors) {
    return MPI_Graph_neighbors_count(comm, rank, nneighbors);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Graphdims_get(AMPI_Comm comm, int* nnodes, int* nedges) {
    return MPI_Graphdims_get(comm, nnodes, nedges);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Topo_test(AMPI_Comm comm, int* status) {
    return MPI_Topo_test(comm, status);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline double AMPI_Wtick() {
    return MPI_Wtick();
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline double AMPI_Wtime() {
    return MPI_Wtime();
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Abort(AMPI_Comm comm, int errorcode) {
    return MPI_Abort(comm, errorcode);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Add_error_class(int* errorclass) {
    return MPI_Add_error_class(errorclass);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Add_error_code(int errorclass, int* errorcode) {
    return MPI_Add_error_code(errorclass, errorcode);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Add_error_string(int errorcode, MEDI_OPTIONAL_CONST char* string) {
    return MPI_Add_error_string(errorcode, string);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Alloc_mem(AMPI_Aint size, AMPI_Info info, void* baseptr) {
    return MPI_Alloc_mem(size, info, baseptr);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_call_errhandler(AMPI_Comm comm, int errorcode) {
    return MPI_Comm_call_errhandler(comm, errorcode);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_create_errhandler(AMPI_Comm_errhandler_function* comm_errhandler_fn, AMPI_Errhandler* errhandler) {
    return MPI_Comm_create_errhandler(comm_errhandler_fn, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_2_2
  inline int AMPI_Comm_create_errhandler(AMPI_Comm_errhandler_fn* comm_errhandler_fn, AMPI_Errhandler* errhandler) {
    return MPI_Comm_create_errhandler(comm_errhandler_fn, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_get_errhandler(AMPI_Comm comm, AMPI_Errhandler* errhandler) {
    return MPI_Comm_get_errhandler(comm, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_set_errhandler(AMPI_Comm comm, AMPI_Errhandler errhandler) {
    return MPI_Comm_set_errhandler(comm, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Errhandler_free(AMPI_Errhandler* errhandler) {
    return MPI_Errhandler_free(errhandler);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Error_class(int errorcode, int* errorclass) {
    return MPI_Error_class(errorcode, errorclass);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Error_string(int errorcode, char* string, int* resultlen) {
    return MPI_Error_string(errorcode, string, resultlen);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_call_errhandler(AMPI_File fh, int errorcode) {
    return MPI_File_call_errhandler(fh, errorcode);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_File_create_errhandler(AMPI_File_errhandler_function* file_errhandler_fn, AMPI_Errhandler* errhandler) {
    return MPI_File_create_errhandler(file_errhandler_fn, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_2_2
  inline int AMPI_File_create_errhandler(AMPI_File_errhandler_fn* file_errhandler_fn, AMPI_Errhandler* errhandler) {
    return MPI_File_create_errhandler(file_errhandler_fn, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_errhandler(AMPI_File file, AMPI_Errhandler* errhandler) {
    return MPI_File_get_errhandler(file, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_set_errhandler(AMPI_File file, AMPI_Errhandler errhandler) {
    return MPI_File_set_errhandler(file, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Finalize() {
    return MPI_Finalize();
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Finalized(int* flag) {
    return MPI_Finalized(flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Free_mem(void* base) {
    return MPI_Free_mem(base);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_Get_library_version(char* version, int* resultlen) {
    return MPI_Get_library_version(version, resultlen);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Get_processor_name(char* name, int* resultlen) {
    return MPI_Get_processor_name(name, resultlen);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Get_version(int* version, int* subversion) {
    return MPI_Get_version(version, subversion);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Initialized(int* flag) {
    return MPI_Initialized(flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_call_errhandler(AMPI_Win win, int errorcode) {
    return MPI_Win_call_errhandler(win, errorcode);
  }

#endif
#if MEDI_MPI_VERSION_2_2 <= MEDI_MPI_TARGET
  inline int AMPI_Win_create_errhandler(AMPI_Win_errhandler_function* win_errhandler_fn, AMPI_Errhandler* errhandler) {
    return MPI_Win_create_errhandler(win_errhandler_fn, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_2_2
  inline int AMPI_Win_create_errhandler(AMPI_Win_errhandler_fn* win_errhandler_fn, AMPI_Errhandler* errhandler) {
    return MPI_Win_create_errhandler(win_errhandler_fn, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_get_errhandler(AMPI_Win win, AMPI_Errhandler* errhandler) {
    return MPI_Win_get_errhandler(win, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Win_set_errhandler(AMPI_Win win, AMPI_Errhandler errhandler) {
    return MPI_Win_set_errhandler(win, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_create(AMPI_Info* info) {
    return MPI_Info_create(info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_delete(AMPI_Info info, MEDI_OPTIONAL_CONST char* key) {
    return MPI_Info_delete(info, key);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_dup(AMPI_Info info, AMPI_Info* newinfo) {
    return MPI_Info_dup(info, newinfo);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_free(AMPI_Info* info) {
    return MPI_Info_free(info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_get(AMPI_Info info, MEDI_OPTIONAL_CONST char* key, int valuelen, char* value, int* flag) {
    return MPI_Info_get(info, key, valuelen, value, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_get_nkeys(AMPI_Info info, int* nkeys) {
    return MPI_Info_get_nkeys(info, nkeys);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_get_nthkey(AMPI_Info info, int n, char* key) {
    return MPI_Info_get_nthkey(info, n, key);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_get_valuelen(AMPI_Info info, MEDI_OPTIONAL_CONST char* key, int* valuelen, int* flag) {
    return MPI_Info_get_valuelen(info, key, valuelen, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Info_set(AMPI_Info info, MEDI_OPTIONAL_CONST char* key, MEDI_OPTIONAL_CONST char* value) {
    return MPI_Info_set(info, key, value);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Close_port(MEDI_OPTIONAL_CONST char* port_name) {
    return MPI_Close_port(port_name);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_accept(MEDI_OPTIONAL_CONST char* port_name, AMPI_Info info, int root, AMPI_Comm comm,
                              AMPI_Comm* newcomm) {
    return MPI_Comm_accept(port_name, info, root, comm, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_connect(MEDI_OPTIONAL_CONST char* port_name, AMPI_Info info, int root, AMPI_Comm comm,
                               AMPI_Comm* newcomm) {
    return MPI_Comm_connect(port_name, info, root, comm, newcomm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_disconnect(AMPI_Comm* comm) {
    return MPI_Comm_disconnect(comm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_get_parent(AMPI_Comm* parent) {
    return MPI_Comm_get_parent(parent);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_join(int fd, AMPI_Comm* intercomm) {
    return MPI_Comm_join(fd, intercomm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_spawn(MEDI_OPTIONAL_CONST char* command, char** argv, int maxprocs, AMPI_Info info, int root,
                             AMPI_Comm comm, AMPI_Comm* intercomm, int* array_of_errcodes) {
    return MPI_Comm_spawn(command, argv, maxprocs, info, root, comm, intercomm, array_of_errcodes);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Comm_spawn_multiple(int count, char** array_of_commands, char*** array_of_argv,
                                      MEDI_OPTIONAL_CONST int* array_of_maxprocs, MEDI_OPTIONAL_CONST AMPI_Info* array_of_info, int root, AMPI_Comm comm,
                                      AMPI_Comm* intercomm, int* array_of_errcodes) {
    return MPI_Comm_spawn_multiple(count, array_of_commands, array_of_argv, array_of_maxprocs, array_of_info, root, comm,
                                   intercomm, array_of_errcodes);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Lookup_name(MEDI_OPTIONAL_CONST char* service_name, AMPI_Info info, char* port_name) {
    return MPI_Lookup_name(service_name, info, port_name);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Open_port(AMPI_Info info, char* port_name) {
    return MPI_Open_port(info, port_name);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Publish_name(MEDI_OPTIONAL_CONST char* service_name, AMPI_Info info,
                               MEDI_OPTIONAL_CONST char* port_name) {
    return MPI_Publish_name(service_name, info, port_name);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Unpublish_name(MEDI_OPTIONAL_CONST char* service_name, AMPI_Info info,
                                 MEDI_OPTIONAL_CONST char* port_name) {
    return MPI_Unpublish_name(service_name, info, port_name);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Grequest_complete(AMPI_Request request) {
    return MPI_Grequest_complete(request.request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Grequest_start(AMPI_Grequest_query_function* query_fn, AMPI_Grequest_free_function* free_fn,
                                 AMPI_Grequest_cancel_function* cancel_fn, void* extra_state, AMPI_Request* request) {
    return MPI_Grequest_start(query_fn, free_fn, cancel_fn, extra_state, &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Is_thread_main(int* flag) {
    return MPI_Is_thread_main(flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Query_thread(int* provided) {
    return MPI_Query_thread(provided);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Status_set_cancelled(AMPI_Status* status, int flag) {
    return MPI_Status_set_cancelled(status, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Status_set_elements(AMPI_Status* status, DATATYPE* datatype, int count) {
    return MPI_Status_set_elements(status, datatype->getModifiedMpiType(), count);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_Status_set_elements_x(AMPI_Status* status, DATATYPE* datatype, AMPI_Count count) {
    return MPI_Status_set_elements_x(status, datatype->getModifiedMpiType(), count);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_close(AMPI_File* fh) {
    return MPI_File_close(fh);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_delete(MEDI_OPTIONAL_CONST char* filename, AMPI_Info info) {
    return MPI_File_delete(filename, info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_amode(AMPI_File fh, int* amode) {
    return MPI_File_get_amode(fh, amode);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_atomicity(AMPI_File fh, int* flag) {
    return MPI_File_get_atomicity(fh, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_byte_offset(AMPI_File fh, AMPI_Offset offset, AMPI_Offset* disp) {
    return MPI_File_get_byte_offset(fh, offset, disp);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_group(AMPI_File fh, AMPI_Group* group) {
    return MPI_File_get_group(fh, group);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_info(AMPI_File fh, AMPI_Info* info_used) {
    return MPI_File_get_info(fh, info_used);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_position(AMPI_File fh, AMPI_Offset* offset) {
    return MPI_File_get_position(fh, offset);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_position_shared(AMPI_File fh, AMPI_Offset* offset) {
    return MPI_File_get_position_shared(fh, offset);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_get_size(AMPI_File fh, AMPI_Offset* size) {
    return MPI_File_get_size(fh, size);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_get_type_extent(AMPI_File fh, DATATYPE* datatype, AMPI_Aint* extent) {
    return MPI_File_get_type_extent(fh, datatype->getModifiedMpiType(), extent);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename ETYPE, typename FILETYPE>
  inline int AMPI_File_get_view(AMPI_File fh, AMPI_Offset* disp, ETYPE* etype, FILETYPE* filetype, char* datarep) {
    return MPI_File_get_view(fh, disp, etype->getModifiedMpiType(), filetype->getModifiedMpiType(), datarep);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iread(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Request* request) {
    return MPI_File_iread(fh, buf, count, datatype->getModifiedMpiType(), request);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iread_all(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Request* request) {
    return MPI_File_iread_all(fh, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iread_at(AMPI_File fh, AMPI_Offset offset, void* buf, int count, DATATYPE* datatype,
                                AMPI_Request* request) {
    return MPI_File_iread_at(fh, offset, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iread_at_all(AMPI_File fh, AMPI_Offset offset, void* buf, int count, DATATYPE* datatype,
                                    AMPI_Request* request) {
    return MPI_File_iread_at_all(fh, offset, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iread_shared(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Request* request) {
    return MPI_File_iread_shared(fh, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iwrite(AMPI_File fh, const void* buf, int count, DATATYPE* datatype, AMPI_Request* request) {
    return MPI_File_iwrite(fh, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iwrite_all(AMPI_File fh, const void* buf, int count, DATATYPE* datatype, AMPI_Request* request) {
    return MPI_File_iwrite_all(fh, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iwrite_at(AMPI_File fh, AMPI_Offset offset, const void* buf, int count, DATATYPE* datatype,
                                 AMPI_Request* request) {
    return MPI_File_iwrite_at(fh, offset, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iwrite_at_all(AMPI_File fh, AMPI_Offset offset, const void* buf, int count, DATATYPE* datatype,
                                     AMPI_Request* request) {
    return MPI_File_iwrite_at_all(fh, offset, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_iwrite_shared(AMPI_File fh, const void* buf, int count, DATATYPE* datatype,
                                     AMPI_Request* request) {
    return MPI_File_iwrite_shared(fh, buf, count, datatype->getModifiedMpiType(), &request->request);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_open(AMPI_Comm comm, MEDI_OPTIONAL_CONST char* filename, int amode, AMPI_Info info,
                            AMPI_File* fh) {
    return MPI_File_open(comm, filename, amode, info, fh);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_preallocate(AMPI_File fh, AMPI_Offset size) {
    return MPI_File_preallocate(fh, size);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_read(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_all(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_read_all(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_all_begin(AMPI_File fh, void* buf, int count, DATATYPE* datatype) {
    return MPI_File_read_all_begin(fh, buf, count, datatype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_read_all_end(AMPI_File fh, void* buf, AMPI_Status* status) {
    return MPI_File_read_all_end(fh, buf, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_at(AMPI_File fh, AMPI_Offset offset, void* buf, int count, DATATYPE* datatype,
                               AMPI_Status* status) {
    return MPI_File_read_at(fh, offset, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_at_all(AMPI_File fh, AMPI_Offset offset, void* buf, int count, DATATYPE* datatype,
                                   AMPI_Status* status) {
    return MPI_File_read_at_all(fh, offset, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_at_all_begin(AMPI_File fh, AMPI_Offset offset, void* buf, int count, DATATYPE* datatype) {
    return MPI_File_read_at_all_begin(fh, offset, buf, count, datatype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_read_at_all_end(AMPI_File fh, void* buf, AMPI_Status* status) {
    return MPI_File_read_at_all_end(fh, buf, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_ordered(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_read_ordered(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_ordered_begin(AMPI_File fh, void* buf, int count, DATATYPE* datatype) {
    return MPI_File_read_ordered_begin(fh, buf, count, datatype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_read_ordered_end(AMPI_File fh, void* buf, AMPI_Status* status) {
    return MPI_File_read_ordered_end(fh, buf, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_read_shared(AMPI_File fh, void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_read_shared(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_seek(AMPI_File fh, AMPI_Offset offset, int whence) {
    return MPI_File_seek(fh, offset, whence);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_seek_shared(AMPI_File fh, AMPI_Offset offset, int whence) {
    return MPI_File_seek_shared(fh, offset, whence);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_set_atomicity(AMPI_File fh, int flag) {
    return MPI_File_set_atomicity(fh, flag);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_set_info(AMPI_File fh, AMPI_Info info) {
    return MPI_File_set_info(fh, info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_set_size(AMPI_File fh, AMPI_Offset size) {
    return MPI_File_set_size(fh, size);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename ETYPE, typename FILETYPE>
  inline int AMPI_File_set_view(AMPI_File fh, AMPI_Offset disp, ETYPE* etype, FILETYPE* filetype, const char* datarep,
                                AMPI_Info info) {
    return MPI_File_set_view(fh, disp, etype->getModifiedMpiType(), filetype->getModifiedMpiType(), datarep, info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_sync(AMPI_File fh) {
    return MPI_File_sync(fh);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write(AMPI_File fh, const void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_write(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_all(AMPI_File fh, const void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_write_all(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_all_begin(AMPI_File fh, const void* buf, int count, DATATYPE* datatype) {
    return MPI_File_write_all_begin(fh, buf, count, datatype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_write_all_end(AMPI_File fh, MEDI_OPTIONAL_CONST void* buf, AMPI_Status* status) {
    return MPI_File_write_all_end(fh, buf, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_at(AMPI_File fh, AMPI_Offset offset, const void* buf, int count, DATATYPE* datatype,
                                AMPI_Status* status) {
    return MPI_File_write_at(fh, offset, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_at_all(AMPI_File fh, AMPI_Offset offset, const void* buf, int count, DATATYPE* datatype,
                                    AMPI_Status* status) {
    return MPI_File_write_at_all(fh, offset, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_at_all_begin(AMPI_File fh, AMPI_Offset offset, const void* buf, int count,
                                          DATATYPE* datatype) {
    return MPI_File_write_at_all_begin(fh, offset, buf, count, datatype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_write_at_all_end(AMPI_File fh, MEDI_OPTIONAL_CONST void* buf, AMPI_Status* status) {
    return MPI_File_write_at_all_end(fh, buf, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_ordered(AMPI_File fh, const void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_write_ordered(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_ordered_begin(AMPI_File fh, const void* buf, int count, DATATYPE* datatype) {
    return MPI_File_write_ordered_begin(fh, buf, count, datatype->getModifiedMpiType());
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_File_write_ordered_end(AMPI_File fh, MEDI_OPTIONAL_CONST void* buf, AMPI_Status* status) {
    return MPI_File_write_ordered_end(fh, buf, status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  template<typename DATATYPE>
  inline int AMPI_File_write_shared(AMPI_File fh, const void* buf, int count, DATATYPE* datatype, AMPI_Status* status) {
    return MPI_File_write_shared(fh, buf, count, datatype->getModifiedMpiType(), status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Register_datarep(MEDI_OPTIONAL_CONST char* datarep,
                                   AMPI_Datarep_conversion_function* read_conversion_fn, AMPI_Datarep_conversion_function* write_conversion_fn,
                                   AMPI_Datarep_extent_function* dtype_file_extent_fn, void* extra_state) {
    return MPI_Register_datarep(datarep, read_conversion_fn, write_conversion_fn, dtype_file_extent_fn, extra_state);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Fint AMPI_Comm_c2f(AMPI_Comm comm) {
    return MPI_Comm_c2f(comm);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Comm AMPI_Comm_f2c(AMPI_Fint comm) {
    return MPI_Comm_f2c(comm);
  }

#endif
#if MEDI_MPI_VERSION_2_1 <= MEDI_MPI_TARGET
  inline MPI_Fint AMPI_Errhandler_c2f(AMPI_Errhandler errhandler) {
    return MPI_Errhandler_c2f(errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_1 <= MEDI_MPI_TARGET
  inline MPI_Errhandler AMPI_Errhandler_f2c(AMPI_Fint errhandler) {
    return MPI_Errhandler_f2c(errhandler);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Fint AMPI_File_c2f(AMPI_File file) {
    return MPI_File_c2f(file);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_File AMPI_File_f2c(AMPI_Fint file) {
    return MPI_File_f2c(file);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Fint AMPI_Group_c2f(AMPI_Group group) {
    return MPI_Group_c2f(group);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Group AMPI_Group_f2c(AMPI_Fint group) {
    return MPI_Group_f2c(group);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Fint AMPI_Info_c2f(AMPI_Info info) {
    return MPI_Info_c2f(info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Info AMPI_Info_f2c(AMPI_Fint info) {
    return MPI_Info_f2c(info);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Status_c2f(MEDI_OPTIONAL_CONST AMPI_Status* c_status, AMPI_Fint* f_status) {
    return MPI_Status_c2f(c_status, f_status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline int AMPI_Status_f2c(MEDI_OPTIONAL_CONST AMPI_Fint* f_status, AMPI_Status* c_status) {
    return MPI_Status_f2c(f_status, c_status);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Fint AMPI_Win_c2f(AMPI_Win win) {
    return MPI_Win_c2f(win);
  }

#endif
#if MEDI_MPI_VERSION_2_0 <= MEDI_MPI_TARGET
  inline MPI_Win AMPI_Win_f2c(AMPI_Fint win) {
    return MPI_Win_f2c(win);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_changed(int* stamp) {
    return MPI_T_category_changed(stamp);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_get_categories(int cat_index, int len, int* indices) {
    return MPI_T_category_get_categories(cat_index, len, indices);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_get_cvars(int cat_index, int len, int* indices) {
    return MPI_T_category_get_cvars(cat_index, len, indices);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_get_index(const char* name, int* cat_index) {
    return MPI_T_category_get_index(name, cat_index);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_get_info(int cat_index, char* name, int* name_len, char* desc, int* desc_len, int* num_cvars,
                                      int* num_pvars, int* num_categories) {
    return MPI_T_category_get_info(cat_index, name, name_len, desc, desc_len, num_cvars, num_pvars, num_categories);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_get_num(int* num_cat) {
    return MPI_T_category_get_num(num_cat);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_category_get_pvars(int cat_index, int len, int* indices) {
    return MPI_T_category_get_pvars(cat_index, len, indices);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  inline int AMPI_T_cvar_get_index(const char* name, int* cvar_index) {
    return MPI_T_cvar_get_index(name, cvar_index);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_cvar_get_num(int* num_cvar) {
    return MPI_T_cvar_get_num(num_cvar);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_cvar_handle_alloc(int cvar_index, void* obj_handle, AMPI_T_cvar_handle* handle, int* count) {
    return MPI_T_cvar_handle_alloc(cvar_index, obj_handle, handle, count);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_cvar_handle_free(AMPI_T_cvar_handle* handle) {
    return MPI_T_cvar_handle_free(handle);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_cvar_read(AMPI_T_cvar_handle handle, void* buf) {
    return MPI_T_cvar_read(handle, buf);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_cvar_write(AMPI_T_cvar_handle handle, const void* buf) {
    return MPI_T_cvar_write(handle, buf);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_enum_get_info(AMPI_T_enum enumtype, int* num, char* name, int* name_len) {
    return MPI_T_enum_get_info(enumtype, num, name, name_len);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_enum_get_item(AMPI_T_enum enumtype, int index, int* value, char* name, int* name_len) {
    return MPI_T_enum_get_item(enumtype, index, value, name, name_len);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_finalize() {
    return MPI_T_finalize();
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_init_thread(int required, int* provided) {
    return MPI_T_init_thread(required, provided);
  }

#endif
#if MEDI_MPI_VERSION_3_1 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_get_index(const char* name, int var_class, int* pvar_index) {
    return MPI_T_pvar_get_index(name, var_class, pvar_index);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_get_num(int* num_pvar) {
    return MPI_T_pvar_get_num(num_pvar);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_handle_alloc(AMPI_T_pvar_session session, int pvar_index, void* obj_handle,
                                      AMPI_T_pvar_handle* handle, int* count) {
    return MPI_T_pvar_handle_alloc(session, pvar_index, obj_handle, handle, count);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_handle_free(AMPI_T_pvar_session session, AMPI_T_pvar_handle* handle) {
    return MPI_T_pvar_handle_free(session, handle);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_read(AMPI_T_pvar_session session, AMPI_T_pvar_handle handle, void* buf) {
    return MPI_T_pvar_read(session, handle, buf);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_readreset(AMPI_T_pvar_session session, AMPI_T_pvar_handle handle, void* buf) {
    return MPI_T_pvar_readreset(session, handle, buf);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_reset(AMPI_T_pvar_session session, AMPI_T_pvar_handle handle) {
    return MPI_T_pvar_reset(session, handle);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_session_create(AMPI_T_pvar_session* session) {
    return MPI_T_pvar_session_create(session);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_session_free(AMPI_T_pvar_session* session) {
    return MPI_T_pvar_session_free(session);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_start(AMPI_T_pvar_session session, AMPI_T_pvar_handle handle) {
    return MPI_T_pvar_start(session, handle);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_stop(AMPI_T_pvar_session session, AMPI_T_pvar_handle handle) {
    return MPI_T_pvar_stop(session, handle);
  }

#endif
#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  inline int AMPI_T_pvar_write(AMPI_T_pvar_session session, AMPI_T_pvar_handle handle, const void* buf) {
    return MPI_T_pvar_write(session, handle, buf);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  inline int AMPI_Address(void* location, AMPI_Aint* address) {
    return MPI_Address(location, address);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Attr_delete(AMPI_Comm comm, int keyval) {
    return MPI_Attr_delete(comm, keyval);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Attr_get(AMPI_Comm comm, int keyval, void* attribute_val, int* flag) {
    return MPI_Attr_get(comm, keyval, attribute_val, flag);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Attr_put(AMPI_Comm comm, int keyval, void* attribute_val) {
    return MPI_Attr_put(comm, keyval, attribute_val);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Keyval_create(AMPI_Copy_function* copy_fn, AMPI_Delete_function* delete_fn, int* keyval,
                                void* extra_state) {
    return MPI_Keyval_create(copy_fn, delete_fn, keyval, extra_state);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
  inline int AMPI_Keyval_free(int* keyval) {
    return MPI_Keyval_free(keyval);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  inline int AMPI_Errhandler_create(AMPI_Handler_function* function, AMPI_Errhandler* errhandler) {
    return MPI_Errhandler_create(function, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  inline int AMPI_Errhandler_get(AMPI_Comm comm, AMPI_Errhandler* errhandler) {
    return MPI_Errhandler_get(comm, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  inline int AMPI_Errhandler_set(AMPI_Comm comm, AMPI_Errhandler errhandler) {
    return MPI_Errhandler_set(comm, errhandler);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  template<typename DATATYPE>
  inline int AMPI_Type_extent(DATATYPE* datatype, AMPI_Aint* extent) {
    return MPI_Type_extent(datatype->getModifiedMpiType(), extent);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  template<typename DATATYPE>
  inline int AMPI_Type_lb(DATATYPE* datatype, AMPI_Aint* displacement) {
    return MPI_Type_lb(datatype->getModifiedMpiType(), displacement);
  }

#endif
#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET && MEDI_MPI_TARGET < MEDI_MPI_VERSION_3_0
  template<typename DATATYPE>
  inline int AMPI_Type_ub(DATATYPE* datatype, AMPI_Aint* displacement) {
    return MPI_Type_ub(datatype->getModifiedMpiType(), displacement);
  }

#endif

#if MEDI_MPI_VERSION_1_0 <= MEDI_MPI_TARGET
#define AMPI_Pcontrol MPI_Pcontrol
#endif


}
