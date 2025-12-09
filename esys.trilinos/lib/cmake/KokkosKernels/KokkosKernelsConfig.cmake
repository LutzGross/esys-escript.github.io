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
# CMake variable for use by Trilinos/KokkosKernels clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "KokkosKernels requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/KokkosKernels build
## ---------------------------------------------------------------------------

set(KokkosKernels_CXX_COMPILER "/usr/bin/g++")

set(KokkosKernels_C_COMPILER "/usr/bin/gcc")

set(KokkosKernels_Fortran_COMPILER "")
# Deprecated!
set(KokkosKernels_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/KokkosKernels build
## ---------------------------------------------------------------------------

## Give the build type
set(KokkosKernels_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(KokkosKernels_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(KokkosKernels_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(KokkosKernels_Fortran_FLAGS [[ ]])
# Deprecated
set(KokkosKernels_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(KokkosKernels_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(KokkosKernels_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(KokkosKernels_BUILD_SHARED_LIBS "ON")

set(KokkosKernels_LINKER /usr/bin/ld)
set(KokkosKernels_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(KokkosKernels_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(KokkosKernels_LIBRARIES KokkosKernels::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(KokkosKernels_MPI_LIBRARIES "")
set(KokkosKernels_MPI_LIBRARY_DIRS "")
set(KokkosKernels_MPI_INCLUDE_DIRS "")
set(KokkosKernels_MPI_EXEC "")
set(KokkosKernels_MPI_EXEC_MAX_NUMPROCS "")
set(KokkosKernels_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(KokkosKernels_ENABLE_Kokkos ON)
set(KokkosKernels_ENABLE_quadmath OFF)
set(KokkosKernels_ENABLE_MKL OFF)
set(KokkosKernels_ENABLE_BLAS ON)
set(KokkosKernels_ENABLE_LAPACK ON)
set(KokkosKernels_ENABLE_METIS OFF)
set(KokkosKernels_ENABLE_SuperLU OFF)
set(KokkosKernels_ENABLE_Cholmod OFF)
set(KokkosKernels_ENABLE_CUBLAS OFF)
set(KokkosKernels_ENABLE_CUSPARSE OFF)
set(KokkosKernels_ENABLE_CUSOLVER OFF)
set(KokkosKernels_ENABLE_ROCBLAS OFF)
set(KokkosKernels_ENABLE_ROCSPARSE OFF)
set(KokkosKernels_ENABLE_ROCSOLVER OFF)

# Exported cache variables
set(KokkosKernels_ENABLE_DEBUG "OFF")
set(HAVE_KOKKOSKERNELS_DEBUG "OFF")

# Include configuration of dependent packages
if (NOT TARGET Kokkos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Kokkos/KokkosConfig.cmake")
endif()
if (NOT TARGET BLAS::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/BLAS/BLASConfig.cmake")
endif()
if (NOT TARGET LAPACK::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/LAPACK/LAPACKConfig.cmake")
endif()

# Import KokkosKernels targets
include("${CMAKE_CURRENT_LIST_DIR}/KokkosKernelsTargets.cmake")

# Standard TriBITS-compliant external package variables
set(KokkosKernels_IS_TRIBITS_COMPLIANT TRUE)
set(KokkosKernels_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(KokkosKernels_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(KokkosKernels_EXPORTED_PACKAGE_LIBS_NAMES "kokkoskernels")

foreach(libname IN LISTS KokkosKernels_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE KokkosKernels::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'KokkosKernels::${libname}', or better yet,"
      " 'KokkosKernels::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'KokkosKernels'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'KokkosKernels_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
