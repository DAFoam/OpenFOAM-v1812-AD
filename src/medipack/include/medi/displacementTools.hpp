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


#include "macros.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

  /**
   * @brief Helper structure for the easy creation of a linear displacement with a
   * the same length.
   */
  struct LinearDisplacements {

      /**
       * @brief The array with the displacements.
       */
      int* displs;

      /**
       * @brief The array with the counts
       */
      int* counts;

      /**
       * @brief Create displacemnts where each displacement has the size length.
       *
       * displs[0] = 0;
       * displs[1] = length;
       * displs[2] = 2 * length;
       * ...
       *
       * counts[i] = length;
       *
       * @param[in] commSize  The size of the communication object.
       * @param[in]   length  The length of each displacement.
       */
      inline LinearDisplacements(int commSize, int length) {
        counts = new int[commSize];
        displs = new int[commSize];
        for(int i = 0; i < commSize; ++i) {
          counts[i] = length;
          displs[i] = i * length;
        }
      }

      /**
       * @brief Destroy the structure.
       */
      inline ~LinearDisplacements() {
        delete [] displs;
        delete [] counts;
      }

      /**
       * @brief Helper function for deleting a LinearDisplacements structure.
       *
       * @param[in] d  The pointer to a LinearDisplacements structure that will be deleted.
       */
      static inline void deleteFunc(void* d) {
        LinearDisplacements* data = reinterpret_cast<LinearDisplacements*>(d);

        delete data;
      }
  };

  /**
   * @brief Compute the total size of a message that has a different size on each rank.
   *
   * @param[in] counts  The size of each rank.
   * @param[in]  ranks  The number of the ranks.
   *
   * @return The sum over all counts.
   */
  inline int computeDisplacementsTotalSize(const int* counts, int ranks) {
    int totalSize = 0;
    for(int i = 0; i < ranks; ++i) {
      totalSize += counts[i];
    }

    return totalSize;
  }

  /**
   * @brief Creates the linearized displacements of a message with a different size on each rank.
   *
   * @param[in] counts  The size of each rank.
   * @param[in]  ranks  The number of the ranks.
   *
   * @return A new displacements array that starts at 0 and increases by the counts on each rank.
   */
  inline int* createLinearDisplacements(const int* counts, int ranks) {
    int* displs = new int[ranks];

    displs[0] = 0;
    for(int i = 1; i < ranks; ++i) {
      displs[i] = counts[i - 1] +  displs[i - 1];
    }

    return displs;
  }

  /**
   * @brief Creates the linearized displacements of a message with a different size on each rank.
   *
   * The counts are scaled by the given factor.
   *
   * @param[in] countsOut  The generated counts.
   * @param[in] displsOut  The generated displacements.
   * @param[in]    counts  The size of each rank.
   * @param[in]     ranks  The number of the ranks.
   * @param[in]     scale  The scaling factor of the counts and displacements.
   */
  inline void createLinearDisplacementsAndCount(int* &countsOut, int* &displsOut, const int* counts, int ranks, int scale) {
    displsOut = new int[ranks];
    countsOut = new int[ranks];

    countsOut[0] = counts[0] * scale;
    displsOut[0] = 0;
    for(int i = 1; i < ranks; ++i) {
      countsOut[i] = counts[i] * scale;
      displsOut[i] = countsOut[i - 1] +  displsOut[i - 1];
    }
  }

  /**
   * @brief Creates the counts for a message with a different size on each rank.
   *
   * The counts are computed by the type such that the result can hold all indices, passive values, etc.
   *
   * @param[out] linearCounts  The linearized counts.
   * @param[in]        counts  The size of each rank.
   * @param[in]        displs  The displacement of each rank.
   * @param[in]         ranks  The number of ranks.
   * @param[in]          type  The mpi data type.
   *
   * @tparam Datatype  The type for the datatype.
   */
  template<typename Datatype>
  inline void createLinearIndexCounts(int* &linearCounts, const int* counts, const int* displs, int ranks, Datatype* type) {
    linearCounts = new int[ranks];

    linearCounts[0] = type->computeActiveElements(displs[0] + counts[0]);
    int lastSum = linearCounts[0];
    for(int i = 1; i < ranks; ++i) {
      int curSum = type->computeActiveElements(displs[i] + counts[i]);
      linearCounts[i] =  curSum - lastSum;
      lastSum = curSum;
    }
  }
}
