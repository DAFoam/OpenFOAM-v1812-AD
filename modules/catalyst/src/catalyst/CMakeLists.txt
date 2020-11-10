#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0002 NEW)   # Policy CMP0002 required for for cmake >= 3

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

include(CMakeLists-OpenFOAM.txt)

#------------------------------------------------------------------------------
# Simple sanity checks

if (EXISTS "$ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED COMPONENTS vtkPVPythonCatalyst)
else()
    message(FATAL_ERROR "ParaView_DIR not set")
endif()


#-----------------------------------------------------------------------------
# All messages

message("================")
message("Using OpenFOAM = $ENV{WM_PROJECT_DIR}")
message("Using ParaView = $ENV{ParaView_DIR}")
if (NOT PARAVIEW_USE_MPI)
    message(WARNING "==== Recommended to build using ParaView Catalyst with MPI ====")
endif()
message("install prefix = ${CMAKE_INSTALL_PREFIX}")
message("================")


#-----------------------------------------------------------------------------

include(CMakeLists-Project.txt)

#-----------------------------------------------------------------------------
