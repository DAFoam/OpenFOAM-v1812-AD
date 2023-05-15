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

#include <mpi.h>

#include "../macros.h"
#include "../typeDefinitions.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

  static void noPreAdjointOperation(void* adjoints, void* primals, int count, int dim) { MEDI_UNUSED(adjoints); MEDI_UNUSED(primals); MEDI_UNUSED(count); MEDI_UNUSED(dim); }
  static void noPostAdjointOperation(void* adjoints, void* primals, void* rootPrimals, int count, int dim) { MEDI_UNUSED(adjoints); MEDI_UNUSED(primals); MEDI_UNUSED(rootPrimals); MEDI_UNUSED(count); MEDI_UNUSED(dim); }

  /**
   * @brief Structure for the special handling of the MPI_Op structure.
   *
   * The structure contains additional data for the AD evaluation and indicators what the
   * operator requires.
   */
  struct AMPI_Op {

      /**
       * @brief Indicates if the primal on the sending and receiving side are required by this operator.
       */
      bool requiresPrimal;

      /**
       * @brief Indicates if the primal on the receiving side needs to be send to the sending side.
       */
      bool requiresPrimalSend;

      /**
       * @brief The mpi operator for the unmodified AD types. The AD tool needs to record all operations
       * that are evaluated with this operator.
       *
       * This operator is used when no specialized adjoint operation is available. Then the buffers are gathered and
       * afterwards the reduction is performed locally.
       */
      MPI_Op primalFunction;

      /**
       * @brief The mpi operator for the modified AD types. This are the operations that are evaluated during the mpi
       * transactions.
       *
       * This operator is used when the user has created a specialized operator for the adjoint communication.
       */
      MPI_Op modifiedPrimalFunction;

      /**
       * @brief The operation that is evaluated on each adjoint value before they are send in a message.
       */
      PreAdjointOperation preAdjointOperation;

      /**
       * @brief The operation that is evaluated on each adjoint after the values have been received in a message.
       */
      PostAdjointOperation postAdjointOperation;

      /**
       * @brief Indicates if the user has provided a specialized adjoint handling for the operator.
       */
      bool hasAdjoint;

      /**
       * @brief Default constructor for static initialization.
       *
       * On this constructed type one of the init method needs to be called.
       */
      AMPI_Op() :
        requiresPrimal(false),
        requiresPrimalSend(false),
        primalFunction(MPI_OP_NULL),
        modifiedPrimalFunction(MPI_OP_NULL),
        preAdjointOperation(noPreAdjointOperation),
        postAdjointOperation(noPostAdjointOperation),
        hasAdjoint(false) {}

      /**
       * @brief Creates an operator with a specialized adjoint handling.
       *
       * @param[in]                requiresPrimal  If the adjoint operations require the primals in the sender/receiver side.
       * @param[in]            requiresPrimalSend  If the adjoint operation on the primal sending side requires the primal result
       *                                           of the primal receiving side.
       * @param[in]                primalFunction  The mpi function for the primal evaluation on the unmodified AD types.
       *                                           This function is called if the buffers are reduced locally and the
       *                                           reduction needs to be recorded by the AD tool.
       * @param[in]         primalFunctionCommute  If the primal function commutes.
       * @param[in]        modifiedPrimalFunction  The mpi function for the primal evaluation of the modifed AD types.
       *                                           This function is called during a MPI transaction, so the AD tool should
       *                                           not record anything.
       * @param[in] modifiedPrimalFunctionCommute  If the modified primal function commutes.
       * @param[in]           preAdjointOperation  The operation that is evaluated on the adjoint values before they are
       *                                           send.
       * @param[in]          postAdjointOperation  The operation that is evaluated on the adjoint values after they are
       *                                           received.
       *
       * @return result of MPI_Op_create
       */
      int init(const bool requiresPrimal, const bool requiresPrimalSend, MPI_User_function* primalFunction, int primalFunctionCommute, MPI_User_function* modifiedPrimalFunction, int modifiedPrimalFunctionCommute, const PreAdjointOperation preAdjointOperation, const PostAdjointOperation postAdjointOperation) {
        int result1 = MPI_Op_create(primalFunction, primalFunctionCommute, &this->primalFunction);
        MPI_Op_create(modifiedPrimalFunction, modifiedPrimalFunctionCommute, &this->modifiedPrimalFunction);

        this->requiresPrimal = requiresPrimal;
        this->requiresPrimalSend = requiresPrimalSend;
        this->preAdjointOperation = preAdjointOperation;
        this->postAdjointOperation = postAdjointOperation;
        this->hasAdjoint = true;

        return result1;
      }

      /**
       * @brief Creates an operator that has no special adjoint handling.
       *
       * If AD types are in this operator, the reduce operations are evaluated locally by performing first a gather
       * operation and the then the reduce.
       *
       * The AD tool needs to record the opertions when the user function is evaluated.
       *
       * @param[in] user_fn  The function for the operation.
       * @param[in] commute  If the operations commutes.
       *
       * @return Result of MPI_Op_Create
       */
      int init(MPI_User_function* user_fn, int commute) {
        MPI_Op op;

        int result = MPI_Op_create(user_fn, commute, &op);
        init(op);

        return result;
      }

      /**
       * @brief Creates an operator that has no special adjoint handling.
       *
       * If AD types are in this operator, the reduce operations are evaluated locally by performing first a gather
       * operation and the then the reduce.
       *
       * The AD tool needs to record the opertions when the user function is evaluated.
       *
       * @param[in] op  An existing MPI operator.
       *
       * @return Result of MPI_Op_Create
       */
      void init(MPI_Op op) {
        this->primalFunction = op;
        this->requiresPrimal = false;
        this->requiresPrimalSend = false;
        this->modifiedPrimalFunction = MPI_SUM;
        this->preAdjointOperation = noPreAdjointOperation;
        this->postAdjointOperation = noPostAdjointOperation;
        this->hasAdjoint = false;
      }

      int free() {
        if(this->hasAdjoint) {
          MPI_Op_free(&this->modifiedPrimalFunction);
        }

        return MPI_Op_free(&this->primalFunction);
      }
  };

  inline bool operator ==(AMPI_Op const& a, AMPI_Op const& b) {
    // Just check the equality for the primal function.
    return a.primalFunction == b.primalFunction;
  }

  inline bool operator !=(AMPI_Op const& a, AMPI_Op const& b) {
    return !(a == b);
  }

  extern const AMPI_Op AMPI_OP_NULL;
}
