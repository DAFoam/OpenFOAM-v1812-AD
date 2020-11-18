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

#include "ampiMisc.h"

#include "../../../generated/medi/ampiDefinitions.h"

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {

#if MEDI_MPI_VERSION_3_0 <= MEDI_MPI_TARGET
  /**
   * @brief Stores additional information for a MPI_Message.
   */
  struct AMPI_Message {
      MPI_Message message;

      int tag;
      int src;
      AMPI_Comm comm;
  };

  inline int AMPI_Mprobe(int source, int tag, AMPI_Comm comm, AMPI_Message* message, AMPI_Status* status) {
    message->src = source;
    message->tag = tag;
    message->comm = comm;

    return MPI_Mprobe(source, tag, comm, &message->message, status);
  }

  inline int AMPI_Improbe(int source, int tag, AMPI_Comm comm, int* flag, AMPI_Message* message, AMPI_Status* status) {
    message->src = source;
    message->tag = tag;
    message->comm = comm;

    return MPI_Improbe(source, tag, comm, flag, &message->message, status);
  }
#endif
}
