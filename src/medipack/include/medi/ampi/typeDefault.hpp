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

#include <new>

#include "../macros.h"
#include "typeInterface.hpp"
#include "op.hpp"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

  /**
   * @brief The default implementation of a MPI type that is represented by an AD type.
   *
   * @tparam ADTool This class needs to implement the ADToolInterface and the StaticADToolInterface.
   */
  template<typename ADTool>
  class MpiTypeDefault final
      : public MpiTypeBase<
          MpiTypeDefault<ADTool>,
          typename ADTool::Type,
          typename ADTool::ModifiedType,
          ADTool>
  {
      INTERFACE_DEF(StaticADToolInterface, ADTool, void)

    public:

      typedef typename ADTool::Type Type;
      typedef typename ADTool::ModifiedType ModifiedType;
      typedef typename ADTool::PrimalType PrimalType;
      typedef typename ADTool::AdjointType AdjointType;
      typedef typename ADTool::IndexType IndexType;

      typedef ADTool Tool;

      bool isClone;

      Tool* adTool;

      MpiTypeDefault(Tool* adTool, MPI_Datatype type, MPI_Datatype modType) :
        MpiTypeBase<MpiTypeDefault<ADTool>, Type, ModifiedType, Tool>(type, modType),
        isClone(false),
        adTool(adTool) {}

    private:
      MpiTypeDefault(Tool* adTool, MPI_Datatype type, MPI_Datatype modType, bool isClone) :
        MpiTypeBase<MpiTypeDefault<ADTool>, Type, ModifiedType, Tool>(type, modType),
        isClone(isClone),
        adTool(adTool) {}

    public:

      ~MpiTypeDefault() {
        if(isClone) {
          MPI_Datatype temp;
          if(this->getModifiedMpiType() != this->getMpiType()) {
            temp = this->getModifiedMpiType();
            MPI_Type_free(&temp);
          }
          temp = this->getMpiType();
          MPI_Type_free(&temp);
        }
      }

      const Tool& getADTool() const {
        return *adTool;
      }

      int computeActiveElements(const int count) const {
        return count;
      }

      bool isModifiedBufferRequired() const {
        return adTool->isModifiedBufferRequired();
      }

      inline void copyIntoModifiedBuffer(const Type* buf, size_t bufOffset, ModifiedType* bufMod, size_t bufModOffset, int elements) const {
        if(adTool->isModifiedBufferRequired()) {
          for(int i = 0; i < elements; ++i) {
            ADTool::setIntoModifyBuffer(bufMod[bufModOffset + i], buf[bufOffset + i]);
          }
        }
      }

      inline void copyFromModifiedBuffer(Type* buf, size_t bufOffset, const ModifiedType* bufMod, size_t bufModOffset, int elements) const {
        if(adTool->isModifiedBufferRequired()) {
          for(int i = 0; i < elements; ++i) {
            ADTool::getFromModifyBuffer(bufMod[bufModOffset + i], buf[bufOffset + i]);
          }
        }
      }

      inline void getIndices(const Type* buf, size_t bufOffset, IndexType* indices, size_t bufModOffset, int elements) const {
        int indexOffset = computeActiveElements((int)bufModOffset);

        for(int i = 0; i < elements; ++i) {
          indices[indexOffset + i] = ADTool::getIndex(buf[bufOffset + i]);
        }
      }

      inline void registerValue(Type* buf, size_t bufOffset, IndexType* indices, PrimalType* oldPrimals, size_t bufModOffset, int elements) const {
        int indexOffset = computeActiveElements((int)bufModOffset);

        for(int i = 0; i < elements; ++i) {
          ADTool::registerValue(buf[bufOffset + i], oldPrimals[indexOffset + i], indices[indexOffset + i]);
        }
      }

      inline void clearIndices(Type* buf, size_t bufOffset, int elements) const {
        for(int i = 0; i < elements; ++i) {
          ADTool::clearIndex(buf[bufOffset + i]);
        }
      }

      inline void createIndices(Type* buf, size_t bufOffset, IndexType* indices, size_t bufModOffset, int elements) const {
        int indexOffset = computeActiveElements((int)bufModOffset);

        for(int i = 0; i < elements; ++i) {
          ADTool::createIndex(buf[bufOffset + i], indices[indexOffset + i]);
        }
      }

      inline void getValues(const Type* buf, size_t bufOffset, PrimalType* primals, size_t bufModOffset, int elements) const {
        int primalOffset = computeActiveElements((int)bufModOffset);

        for(int pos = 0; pos < elements; ++pos) {
          primals[primalOffset + pos] = ADTool::getValue(buf[bufOffset + pos]);
        }
      }

      inline void performReduce(Type* buf, Type* target, int count, AMPI_Op op, int ranks) const {
        for(int j = 1; j < ranks; ++j) {
          MPI_Reduce_local(&buf[j * count], buf, count, this->getMpiType(), op.primalFunction);
        }

        if(0 != ranks) {
          copy(buf, 0, target, 0, count);
        }
      }

      void copy(Type* from, size_t fromOffset, Type* to, size_t toOffset, int count) const {
        for(int i = 0; i < count; ++i) {
          to[toOffset + i] = from[fromOffset + i];
        }

      }

      void initializeType(Type* buf, size_t bufOffset, int elements) const {
        for(int i = 0; i < elements; ++i) {
          new(&buf[bufOffset + i]) Type;
        }
      }

      void freeType(Type* buf, size_t bufOffset, int elements) const {
        for(int i = 0; i < elements; ++i) {
          buf[bufOffset + i].~Type();
        }
      }

      inline void createTypeBuffer(Type* &buf, size_t size) const {
        buf = new Type[size];
      }

      inline void createModifiedTypeBuffer(ModifiedType* &buf, size_t size) const {
        buf = new ModifiedType[size];
      }

      inline void deleteTypeBuffer(Type* &buf, size_t size) const {
        MEDI_UNUSED(size);

        if(NULL != buf) {
          delete [] buf;
          buf = NULL;
        }
      }

      inline void deleteModifiedTypeBuffer(ModifiedType* &buf) const {
        if(NULL != buf) {
          delete [] buf;
          buf = NULL;
        }
      }

      inline MpiTypeDefault* clone() const {
        MPI_Datatype type;
        MPI_Datatype modType;

        MPI_Type_dup(this->getMpiType(), &type);
        if(this->getMpiType() != this->getModifiedMpiType()) {
          MPI_Type_dup(this->getModifiedMpiType(), &modType);
        } else {
          modType = type;
        }

        return new MpiTypeDefault(adTool, this->getMpiType(), this->getModifiedMpiType(), true);
      }
  };
}
