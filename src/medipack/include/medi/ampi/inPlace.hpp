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

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

  /**
   * @brief Implementation of the MPI_IN_PLACE structure.
   *
   * This structures implements cast operators such that it can be used by any type.
   */
  struct AMPI_IN_PLACE_IMPL {
    
      AMPI_IN_PLACE_IMPL () {}

      template<typename T>
      operator const T*() const {
          return reinterpret_cast<const T*>(MPI_IN_PLACE);
      }

      template<typename T>
      operator T*() const {
          return reinterpret_cast<T*>(MPI_IN_PLACE);
      }

  };

  /**
   * @brief This structure is able to be cast to any type.
   */
  extern const AMPI_IN_PLACE_IMPL AMPI_IN_PLACE;
}
