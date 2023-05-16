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

#include "adToolInterface.h"
#include "macros.h"
#include "typeDefinitions.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

  /**
   * @brief Implementation for the AD tool interface of a type that is no AD type.
   *
   * All methods in this implementation contain no logic.
   */
  class ADToolPassive final : public ADToolBase<ADToolPassive, void, void, void> {
    public:

      typedef void PrimalType;
      typedef void AdjointType;
      typedef void IndexType;

      ADToolPassive(MPI_Datatype primalType, MPI_Datatype adjointType) :
        ADToolBase<ADToolPassive, void, void, void>(primalType, adjointType)
      {}

      inline bool isActiveType() const {return false;}
      inline bool isHandleRequired() const {return false;}
      inline bool isModifiedBufferRequired() const {return false;}
      inline bool isOldPrimalsRequired() const {return false;}
      inline void startAssembly(HandleBase* h) const {MEDI_UNUSED(h);}
      inline void stopAssembly(HandleBase* h) const {MEDI_UNUSED(h);}
      inline void addToolAction(HandleBase* h) const {MEDI_UNUSED(h);}

      inline AMPI_Op convertOperator(AMPI_Op op) const {
        return op;
      }

      inline void createPrimalTypeBuffer(PrimalType* &buf, size_t size) const {
        MEDI_UNUSED(size);

        buf = nullptr;
      }

      inline void createIndexTypeBuffer(IndexType* &buf, size_t size) const {
        MEDI_UNUSED(size);

        buf = nullptr;
      }

      inline void deletePrimalTypeBuffer(PrimalType* &buf) const {
        buf = nullptr;
      }

      inline void deleteIndexTypeBuffer(IndexType* &buf) const {
        buf = nullptr;
      }
  };
}
