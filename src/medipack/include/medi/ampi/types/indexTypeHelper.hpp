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

#include <algorithm>

#include "../ampiMisc.h"
#include "../../macros.h"
#include "../typeInterface.hpp"
#include "../op.hpp"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

  template <typename Type, typename ModifiedType, typename PrimalType, typename IndexType>
  struct ToolInterface {
      static PrimalType getPrimalFromMod(const ModifiedType& mod);
      static void setPrimalToMod(ModifiedType& mod, const PrimalType& value);

      static void modifyDependency(const ModifiedType& in, ModifiedType& inout);
  };

  /**
   * @brief The provides all methods required for the creation of operators for AD types.
   *
   * @tparam          Type  The floating point type of the AD tool
   * @tparam  ModifiedType  The type that is send over the network.
   * @tparam    PrimalType  The primal floating point type which is replaced by the AD type.
   * @tparam     IndexType  The identifier used by the AD tool for the AD types.
   * @tparam          Tool  The interface to the AD tool required by this class. The type needs to implement the ToolInterface class.
   */
  template <typename Type, typename ModifiedType, typename PrimalType, typename IndexType, typename AdjointType, INTERFACE_ARG(Tool)>
  struct FunctionHelper {
      INTERFACE_DEF(ToolInterface, Tool, Type, ModifiedType, PrimalType, IndexType)

      struct TypeInt {
          Type value;
          int index;
      };

      struct ModTypeInt {
          ModifiedType value;
          int index;
      };

      static void unmodifiedAdd(Type* invec, Type* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        for(int i = 0; i < *len; ++i) {
          inoutvec[i] += invec[i];
        }
      }

      static void unmodifiedMul(Type* invec, Type* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        for(int i = 0; i < *len; ++i) {
          inoutvec[i] *= invec[i];
        }
      }

      static void unmodifiedMax(Type* invec, Type* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::max;
        for(int i = 0; i < *len; ++i) {
          inoutvec[i] = max(inoutvec[i], invec[i]);
        }
      }

      static void unmodifiedMin(Type* invec, Type* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::min;
        for(int i = 0; i < *len; ++i) {
          inoutvec[i] = min(inoutvec[i], invec[i]);
        }
      }

      static void unmodifiedMaxLoc(TypeInt* invec, TypeInt* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::max;
        for(int i = 0; i < *len; ++i) {
          // first determine the index
          if(invec[i].value > inoutvec[i].value) {
            inoutvec[i].index = invec[i].index;
          } else if(invec[i].value < inoutvec[i].value){
            // empty operation: inoutvec[i].index = inoutvec[i].index;
          } else {
            inoutvec[i].index = std::min(invec[i].index, inoutvec[i].index);
          }

          // second determine the value
          inoutvec[i].value = max(inoutvec[i].value, invec[i].value);
        }
      }

      static void unmodifiedMinLoc(TypeInt* invec, TypeInt* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::min;
        for(int i = 0; i < *len; ++i) {
          // first determine the index
          if(invec[i].value < inoutvec[i].value) {
            inoutvec[i].index = invec[i].index;
          } else if(invec[i].value > inoutvec[i].value){
            // empty operation: inoutvec[i].index = inoutvec[i].index;
          } else {
            inoutvec[i].index = std::min(invec[i].index, inoutvec[i].index);
          }

          // second determine the value
          inoutvec[i].value = min(inoutvec[i].value, invec[i].value);
        }
      }

      static void modifiedAdd(ModifiedType* invec, ModifiedType* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        for(int i = 0; i < *len; ++i) {
          Tool::modifyDependency(invec[i], inoutvec[i]);
          Tool::setPrimalToMod(inoutvec[i], Tool::getPrimalFromMod(invec[i]) + Tool::getPrimalFromMod(inoutvec[i]));
        }
      }

      static void modifiedMul(ModifiedType* invec, ModifiedType* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        for(int i = 0; i < *len; ++i) {
          Tool::modifyDependency(invec[i], inoutvec[i]);
          Tool::setPrimalToMod(inoutvec[i], Tool::getPrimalFromMod(invec[i]) * Tool::getPrimalFromMod(inoutvec[i]));
        }
      }

      static void modifiedMax(ModifiedType* invec, ModifiedType* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::max;
        for(int i = 0; i < *len; ++i) {
          Tool::modifyDependency(invec[i], inoutvec[i]);
          Tool::setPrimalToMod(inoutvec[i], max(Tool::getPrimalFromMod(invec[i]), Tool::getPrimalFromMod(inoutvec[i])));
        }
      }

      static void modifiedMin(ModifiedType* invec, ModifiedType* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::min;
        for(int i = 0; i < *len; ++i) {
          Tool::modifyDependency(invec[i], inoutvec[i]);
          Tool::setPrimalToMod(inoutvec[i], min(Tool::getPrimalFromMod(invec[i]), Tool::getPrimalFromMod(inoutvec[i])));
        }
      }

      static void modifiedMaxLoc(ModTypeInt* invec, ModTypeInt* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::max;
        for(int i = 0; i < *len; ++i) {

          PrimalType inPrimal = Tool::getPrimalFromMod(invec[i].value);
          PrimalType inoutPrimal = Tool::getPrimalFromMod(inoutvec[i].value);

          // first determine the index
          if(inPrimal > inoutPrimal) {
            inoutvec[i].index = invec[i].index;
          } else if(inPrimal < inoutPrimal){
            // empty operation: inoutvec[i].index = inoutvec[i].index;
          } else {
            inoutvec[i].index = std::min(invec[i].index, inoutvec[i].index);
          }

          Tool::modifyDependency(invec[i].value, inoutvec[i].value);
          Tool::setPrimalToMod(inoutvec[i].value, max(inPrimal, inoutPrimal));
        }
      }

      static void modifiedMinLoc(ModTypeInt* invec, ModTypeInt* inoutvec, int* len, MPI_Datatype* datatype) {
        MEDI_UNUSED(datatype);

        using std::min;
        for(int i = 0; i < *len; ++i) {
          PrimalType inPrimal = Tool::getPrimalFromMod(invec[i].value);
          PrimalType inoutPrimal = Tool::getPrimalFromMod(inoutvec[i].value);

          // first determine the index
          if(inPrimal < inoutPrimal) {
            inoutvec[i].index = invec[i].index;
          } else if(inPrimal > inoutPrimal){
            // empty operation: inoutvec[i].index = inoutvec[i].index;
          } else {
            inoutvec[i].index = std::min(invec[i].index, inoutvec[i].index);
          }

          Tool::modifyDependency(invec[i].value, inoutvec[i].value);
          Tool::setPrimalToMod(inoutvec[i].value, min(inPrimal, inoutPrimal));
        }
      }

//      TODO: These are currently not used since we can not handle zero terms. Need to implement
//            a tracking of how many zeros the multiplication contained.
//
//      void preAdjMul(AdjointType* adjoints, PrimalType* primals, int count) {
//        for(int i = 0; i < count; ++i) {
//          adjoints[i] *= primals[i];
//        }
//      }
//
//      void postAdjMul(AdjointType* adjoints, PrimalType* primals, PrimalType* rootPrimals, int count) {
//        CODI_UNUSED(rootPrimals);
//
//        for(int i = 0; i < count; ++i) {
//          if(0.0 != primals[i]) {
//            adjoints[i] /= primals[i];
//          }
//        }
//      }

      static void postAdjMinMax(AdjointType* adjoints, PrimalType* primals, PrimalType* rootPrimals, int count, int vecSize) {
        for(int i = 0; i < count; ++i) {
          if(rootPrimals[i] != primals[i]) {
            for(int dim = 0; dim < vecSize; ++dim) {
              adjoints[i * vecSize + dim] = AdjointType(); // the primal of this process was not the minimum or maximum so do not perfrom the adjoint update
            }
          }
        }
      }
  };

  /**
   *
   * @tparam ADTypeImpl  Needs to implement the adToolInferface and provide the static elements for the type.
   */
  template<INTERFACE_ARG(FuncHelp)>
  struct OperatorHelper {
    public:

      INTERFACE_DEF(FunctionHelper, FuncHelp, void, void, void, void)

      AMPI_Op OP_SUM;
      AMPI_Op OP_PROD;
      AMPI_Op OP_MIN;
      AMPI_Op OP_MAX;
      AMPI_Op OP_MINLOC;
      AMPI_Op OP_MAXLOC;

      void createOperators() {
        AMPI_Op_create((MPI_User_function*)FuncHelp::unmodifiedAdd, 1, &OP_SUM);
        AMPI_Op_create((MPI_User_function*)FuncHelp::unmodifiedMul, 1, &OP_PROD);
        AMPI_Op_create((MPI_User_function*)FuncHelp::unmodifiedMin, 1, &OP_MIN);
        AMPI_Op_create((MPI_User_function*)FuncHelp::unmodifiedMax, 1, &OP_MAX);
        AMPI_Op_create((MPI_User_function*)FuncHelp::unmodifiedMinLoc, 1, &OP_MINLOC);
        AMPI_Op_create((MPI_User_function*)FuncHelp::unmodifiedMaxLoc, 1, &OP_MAXLOC);
      }

      AMPI_Op convertOperator(AMPI_Op op) const {
        if(MPI_SUM == op.primalFunction) {
          return OP_SUM;
        } else if(MPI_PROD == op.primalFunction) {
          return OP_PROD;
        } else if(MPI_MIN == op.primalFunction) {
          return OP_MIN;
        } else if(MPI_MAX == op.primalFunction) {
          return OP_MAX;
        } else if(MPI_MINLOC == op.primalFunction) {
          return OP_MINLOC;
        } else if(MPI_MAXLOC == op.primalFunction) {
          return OP_MAXLOC;
        } else {
          // do not change the type if it is not one of the above
          return op;
        }
      }

      void init() {
        createOperators();
      }

      void finalize() {
        OP_SUM.free();
        OP_PROD.free();
        OP_MIN.free();
        OP_MAX.free();
        OP_MINLOC.free();
        OP_MAXLOC.free();
      }

      static AMPI_Datatype createIntType(const AMPI_Datatype type) {

        AMPI_Datatype intType;
        AMPI_Aint offsets[3] = {
          offsetof(typename FuncHelp::TypeInt, value),
          offsetof(typename FuncHelp::TypeInt, index),
          offsetof(typename FuncHelp::TypeInt, index) + sizeof(int)
        };
        int blockLength[3] = {1, 1, (int)(sizeof(typename FuncHelp::TypeInt) - offsets[2])};
        const AMPI_Datatype types[3] = {type, AMPI_INT, AMPI_BYTE};

        AMPI_Type_create_struct(3, blockLength, offsets, types, &intType);
        AMPI_Type_commit(&intType);

        return intType;
      }

      static void freeIntType(AMPI_Datatype& type) {
        AMPI_Type_free(&type);
      }
  };
}
