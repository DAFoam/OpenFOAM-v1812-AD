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

#include "typeDefinitions.h"
#include "adToolInterface.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {


  class AdjointInterface {
    public:

      /**
       * @brief Compute the number of active types in the buffer.
       * @param[in] elements  The number of elements in the buffer.
       * @return The number of active types in the buffer.
       */
      virtual int computeElements(int elements) const = 0;

      /**
       * @brief The vector size for the current evaluation.
       * @return The vector size for the current evaluation.
       */
      virtual int getVectorSize() const = 0;

      /**
       * @brief Create an array for the primal variables.
       *
       * @param[out] buf  The pointer for the buffer.
       * @param[in] size  The size of the buffer.
       */
      virtual void createPrimalTypeBuffer(void* &buf, size_t size) const = 0;

      /**
       * @brief Delete the array of the primal variables.
       *
       * @param[in,out] buf  The pointer for the buffer.
       */
      virtual void deletePrimalTypeBuffer(void* &buf) const = 0;

      /**
       * @brief Create an array for the adjoint variables.
       *
       * @param[out] buf  The pointer for the buffer.
       * @param[in] size  The size of the buffer.
       */
      virtual void createAdjointTypeBuffer(void* &buf, size_t size) const = 0;

      /**
       * @brief Delete the array of the adjoint variables.
       *
       * @param[in,out] buf  The pointer for the buffer.
       */
      virtual void deleteAdjointTypeBuffer(void* &buf) const = 0;

      /**
       * @brief Perform a reduction in the first element of the buffer.
       * @param[in,out]  buf  The buffer with adjoint values its size is elements * ranks
       * @param[in] elements  The number of elements in the vectors.
       * @param[in]    ranks  The number of ranks in the communication.
       */
      virtual void combineAdjoints(void* buf, const int elements, const int ranks) const = 0;

      /**
       * @brief Get the adjoints for the indices from the AD tool.
       *
       * @param[in]   indices  The indices from the AD tool for the variables in the buffer.
       * @param[out] adjoints  The vector for the adjoint variables.
       * @param[in]  elements  The number of elements in the vectors.
       */
      virtual void getAdjoints(const void* indices, void* adjoints, int elements) const = 0;

      /**
       * @brief Add the adjoint varaibles to the ones in the AD tool. That is the AD tool should perform the
       * operation:
       *
       * internalAdjoints[indices[i]] += adjoints[i];
       *
       * @param[in]   indices  The indices from the AD tool for the variables in the buffer.
       * @param[out] adjoints  The vector with the adjoint variables.
       * @param[in]  elements  The number of elements in the vectors.
       */
      virtual void updateAdjoints(const void* indices, const void* adjoints, int elements) const = 0;

      /**
       * @brief Get the primal values from the AD tool.
       *
       * Can be used to store the old primal values from the floating point values in the buffer.
       *
       * @param[in]   indices  The indices from the AD tool for the variables in the buffer.
       * @param[out]  primals  The vector with the old primal variables.
       * @param[in]  elements  The number of elements in the vectors.
       */
      virtual void getPrimals(const void* indices, const void* primals, int elements) const = 0;

      /**
       * @brief Set the primal values on the AD tool.
       *
       * Can be used to restore the old primal values from the floating point values in the buffer.
       *
       * @param[in]   indices  The indices from the AD tool for the variables in the buffer.
       * @param[out]  primals  The vector with the old primal variables.
       * @param[in]  elements  The number of elements in the vectors.
       */
      virtual void setPrimals(const void* indices, const void* primals, int elements) const = 0;
  };
}
