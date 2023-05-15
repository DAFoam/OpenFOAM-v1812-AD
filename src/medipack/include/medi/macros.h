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

#define MEDI_MPI_VERSION_1_0 100
#define MEDI_MPI_VERSION_1_1 101
#define MEDI_MPI_VERSION_1_2 102
#define MEDI_MPI_VERSION_1_3 103
#define MEDI_MPI_VERSION_2_0 200
#define MEDI_MPI_VERSION_2_1 201
#define MEDI_MPI_VERSION_2_2 202
#define MEDI_MPI_VERSION_3_0 300
#define MEDI_MPI_VERSION_3_1 301


#ifndef MEDI_MPI_TARGET
# if defined(MPI_VERSION) && defined(MPI_SUBVERSION)
#   define MEDI_MPI_TARGET (MPI_VERSION * 100 + MPI_SUBVERSION)
# else
#   warning MEDI: Could not detect MPI version please define MPI_VERSION and MPI_SUBVERSION. Using MPI 3.1 as default.
#   define MEDI_MPI_TARGET MEDI_MPI_VERSION_3_1
# endif
#endif

#if !defined MEDI_NO_OPTIONAL_CONST
# if defined (OPEN_MPI)
#   define OMPI_VERSION (OMPI_MAJOR_VERSION * 100 + OMPI_MINOR_VERSION * 10)
#   if (OMPI_VERSION < 170)
#     define MEDI_OPTIONAL_CONST /* const */
#   else
#     define MEDI_OPTIONAL_CONST const
#   endif
# elif defined (MPICH2)
#   if (MPICH2_NUMVERSION < 15000000)
#     define MEDI_OPTIONAL_CONST /* const */
#   elif (MPICH2_NUMVERSION < 30000000)
#     define MEDI_OPTIONAL_CONST MPICH2_CONST
#   else
#     define MEDI_OPTIONAL_CONST const
#   endif
# else
#   define MEDI_OPTIONAL_CONST const
# endif
#else
# define MEDI_OPTIONAL_CONST /* const */
#endif


#ifndef MEDI_EnableAssert
  #define MEDI_EnableAssert 1
#endif
#ifndef mediAssert
  #if MEDI_EnableAssert

#include <assert.h>
    /**
     * @brief The assert function for MeDiPack it can be enabled with the preprocessor macro MEDI_EnableAssert=true
     *
     * @param x The expression that is checked in the assert.
     *
     * It can be set with the preprocessor macro MEDI_EnableAssert=<true/false>
     */
    #define mediAssert(x) assert(x)
  #else
    /**
     * @brief The assert function for MeDiPack it can be enabled with the preprocessor macro MEDI_EnableAssert=true
     *
     * It can be set with the preprocessor macro MEDI_EnableAssert=<true/false>
     *
     * @param x The expression that is checked in the assert.
     */
    #define mediAssert(x) /* disabled by MEDI_EnableAssert */
  #endif
#endif

/**
 * @brief Global namespace for MeDiPack - Message Differentiation Package
 */
namespace medi {
  #define MEDI_UNUSED(name) (void)(name)
  #define MEDI_CHECK_ERROR(expr) (expr)

  #ifdef DEV
    #define INTERFACE_ARG(name) bool __p
    #define INTERFACE_DEF(interface, name, ...) typedef interface<__VA_ARGS__> name;
  #else
    #define INTERFACE_ARG(name) typename name
    #define INTERFACE_DEF(interface, name, ...) /* empty */
  #endif
}
