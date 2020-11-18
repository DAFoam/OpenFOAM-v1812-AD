/*
 * MeDiPack, a Message Differentiation Package
 *
 * Copyright (C) 2017-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "../macros.h"
#include "op.hpp"
#include "../../../generated/medi/ampiDefinitions.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {
  /**
   * @brief Default forward of the operator creation.
   *
   * Operators created with this function will not be able to perform optimized reduction operations.
   * See the overloads for the creation of operators, that have optimized reduction.
   *
   * For this operator reduce operations are changed into gather operations and then the reduce is computed locally.
   *
   * @param[in] user_fn  The user function. This function needs to use the AD types for the operations.
   * @param[in] commute  If the operation commutes.
   * @param[out]     op  The address of the operator.
   *
   * @return MPI status code.
   */
  inline int AMPI_Op_create(MPI_User_function* user_fn, int commute, AMPI_Op* op) {
    return op->init(user_fn, commute);
  }

  /**
   * @brief Optimized creation of an operator that will enalbe the AD tool to execute reduce operations.
   *
   * The primal function is the same as the usual function and has to use the AD type.
   * The modified primal function performs the same operation as the primal function but with the
   * modified type that the AD tool uses in the buffers.
   *
   * The pre adjoint operations is executed prior to sending the adjoint values to the original sender side.
   * The post adjoint operation is executed bevor the adjoitn values are updated on the sender side.
   *
   * TODO: Add latex description.
   *
   * @param[in]                requiresPrimal  If the primal values should be stored and made available to the pre and post operations. Otherwise the pointers in thes function calls will be zero.
   * @param[in]            requiresPrimalSend  If the sender side requires the primal values of the reciever side for the pre and post operations. Otherwise the pointers in thes function calls will be zero.
   * @param[in]                primalFunction  The function that evaluates the operation on the AD types.
   * @param[in]         primalFunctionCommute  If the function commutes.
   * @param[in]        modifiedPrimalFunction  The function that evaluates the operation on the modified AD types.
   * @param[in] modifiedPrimalFunctionCommute  If the function commutes.
   * @param[in]           preAdjointOperation  The operation that is evaluated before the adjoint values are send to the sender side.
   * @param[in]          postAdjointOperation  The operation that is evaluated before the adjoint values are updated on the sender side.
   * @param[out]                           op  The address of the operator.
   *
   * @return MPI status code.
   */
  inline int AMPI_Op_create(const bool requiresPrimal, const bool requiresPrimalSend,
                     MPI_User_function* primalFunction, int primalFunctionCommute,
                     MPI_User_function* modifiedPrimalFunction, int modifiedPrimalFunctionCommute,
                     const PreAdjointOperation preAdjointOperation,
                     const PostAdjointOperation postAdjointOperation,
                     AMPI_Op* op) {
    return op->init(requiresPrimal, requiresPrimalSend, primalFunction, primalFunctionCommute, modifiedPrimalFunction, modifiedPrimalFunctionCommute, preAdjointOperation, postAdjointOperation);
  }

  /**
   * @brief Frees the operator.
   *
   * @param[in,out] op  The address of the operator.
   *
   * @return  MPI status code.
   */
  inline int AMPI_Op_free(AMPI_Op* op) {
    return op->free();
  }
}
