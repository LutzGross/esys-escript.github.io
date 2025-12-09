# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

##############################################################################
#
# CMake variable for use by Trilinos/Ifpack2 clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Ifpack2 requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Ifpack2 build
## ---------------------------------------------------------------------------

set(Ifpack2_CXX_COMPILER "/usr/bin/g++")

set(Ifpack2_C_COMPILER "/usr/bin/gcc")

set(Ifpack2_Fortran_COMPILER "")
# Deprecated!
set(Ifpack2_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Ifpack2 build
## ---------------------------------------------------------------------------

## Give the build type
set(Ifpack2_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Ifpack2_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Ifpack2_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(Ifpack2_Fortran_FLAGS [[ ]])
# Deprecated
set(Ifpack2_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(Ifpack2_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Ifpack2_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(Ifpack2_BUILD_SHARED_LIBS "ON")

set(Ifpack2_LINKER /usr/bin/ld)
set(Ifpack2_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Ifpack2_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(Ifpack2_LIBRARIES Ifpack2::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Ifpack2_MPI_LIBRARIES "")
set(Ifpack2_MPI_LIBRARY_DIRS "")
set(Ifpack2_MPI_INCLUDE_DIRS "")
set(Ifpack2_MPI_EXEC "")
set(Ifpack2_MPI_EXEC_MAX_NUMPROCS "")
set(Ifpack2_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Ifpack2_ENABLE_Belos ON)
set(Ifpack2_ENABLE_Teuchos ON)
set(Ifpack2_ENABLE_Tpetra ON)
set(Ifpack2_ENABLE_KokkosKernels ON)
set(Ifpack2_ENABLE_Xpetra ON)
set(Ifpack2_ENABLE_Zoltan2Core OFF)
set(Ifpack2_ENABLE_ThyraTpetraAdapters OFF)
set(Ifpack2_ENABLE_Amesos2 ON)
set(Ifpack2_ENABLE_ShyLU_NodeBasker OFF)
set(Ifpack2_ENABLE_ShyLU_NodeHTS OFF)
set(Ifpack2_ENABLE_ShyLU_NodeFastILU OFF)
set(Ifpack2_ENABLE_HYPRE OFF)
set(Ifpack2_ENABLE_Cholmod OFF)
set(Ifpack2_ENABLE_Lemon OFF)
set(Ifpack2_ENABLE_METIS OFF)
set(Ifpack2_ENABLE_MPI OFF)

# Exported cache variables
set(Ifpack2_Trilinos "ON")
set(HAVE_Ifpack2_Trilinos "ON")
set(Ifpack2_ENABLE_DEBUG "OFF")
set(HAVE_IFPACK2_DEBUG "OFF")
set(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_IFPACK2_EXPLICIT_INSTANTIATION "ON")
set(Ifpack2_ENABLE_Experimental "NO")
set(HAVE_IFPACK2_EXPERIMENTAL "OFF")
set(Ifpack2_ENABLE_DEPRECATED_CODE "NO")
set(HAVE_IFPACK2_DEPRECATED_CODE "OFF")
set(Ifpack2_ENABLE_Experimental_KokkosKernels_Features "ON")
set(HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES "ON")
set(Ifpack2_ENABLE_IFPACK2_TIMER_BARRIER "NO")
set(HAVE_IFPACK2_TIMER_BARRIER "OFF")
set(Ifpack2_ENABLE_BlockTriDiContainer_Timers "NO")
set(HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS "OFF")
set(Ifpack2_ENABLE_BlockTriDiContainer_Small_Scalar "OFF")
set(HAVE_IFPACK2_BLOCKTRIDICONTAINER_SMALL_SCALAR "OFF")

# Include configuration of dependent packages
if (NOT TARGET Belos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Belos/BelosConfig.cmake")
endif()
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Tpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Tpetra/TpetraConfig.cmake")
endif()
if (NOT TARGET KokkosKernels::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../KokkosKernels/KokkosKernelsConfig.cmake")
endif()
if (NOT TARGET Xpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Xpetra/XpetraConfig.cmake")
endif()
if (NOT TARGET Amesos2::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Amesos2/Amesos2Config.cmake")
endif()

# Import Ifpack2 targets
include("${CMAKE_CURRENT_LIST_DIR}/Ifpack2Targets.cmake")

# Standard TriBITS-compliant external package variables
set(Ifpack2_IS_TRIBITS_COMPLIANT TRUE)
set(Ifpack2_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Ifpack2_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Ifpack2_EXPORTED_PACKAGE_LIBS_NAMES "ifpack2")

foreach(libname IN LISTS Ifpack2_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Ifpack2::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Ifpack2::${libname}', or better yet,"
      " 'Ifpack2::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Ifpack2'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Ifpack2_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
