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
# CMake variable for use by Trilinos/Tpetra clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Tpetra requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Tpetra build
## ---------------------------------------------------------------------------

set(Tpetra_CXX_COMPILER "/usr/bin/g++")

set(Tpetra_C_COMPILER "/usr/bin/gcc")

set(Tpetra_Fortran_COMPILER "")
# Deprecated!
set(Tpetra_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Tpetra build
## ---------------------------------------------------------------------------

## Give the build type
set(Tpetra_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Tpetra_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Tpetra_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(Tpetra_Fortran_FLAGS [[ ]])
# Deprecated
set(Tpetra_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(Tpetra_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Tpetra_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(Tpetra_BUILD_SHARED_LIBS "ON")

set(Tpetra_LINKER /usr/bin/ld)
set(Tpetra_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Tpetra_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(Tpetra_LIBRARIES Tpetra::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Tpetra_MPI_LIBRARIES "")
set(Tpetra_MPI_LIBRARY_DIRS "")
set(Tpetra_MPI_INCLUDE_DIRS "")
set(Tpetra_MPI_EXEC "")
set(Tpetra_MPI_EXEC_MAX_NUMPROCS "")
set(Tpetra_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Tpetra_ENABLE_TpetraCore ON)
set(Tpetra_ENABLE_TpetraTSQR ON)

# Exported cache variables
set(Tpetra_ENABLE_DEBUG "OFF")
set(HAVE_TPETRA_DEBUG "OFF")
set(Tpetra_ENABLE_DEPRECATED_CODE "ON")
set(TPETRA_ENABLE_DEPRECATED_CODE "ON")
set(Tpetra_IGNORE_KOKKOS_COMPATIBILITY "OFF")
set(Tpetra_ENABLE_BLOCKCRS_LITTLEBLOCK_LAYOUTLEFT "OFF")
set(TPETRA_ENABLE_BLOCKCRS_LITTLEBLOCK_LAYOUTLEFT "OFF")
set(Tpetra_ENABLE_CUDA "OFF")
set(TPETRA_ENABLE_CUDA "OFF")
set(Tpetra_INST_SERIAL "ON")
set(HAVE_TPETRA_INST_SERIAL "ON")
set(Tpetra_INST_OPENMP "ON")
set(HAVE_TPETRA_INST_OPENMP "ON")
set(Tpetra_INST_PTHREAD "OFF")
set(HAVE_TPETRA_INST_PTHREAD "OFF")
set(Tpetra_INST_CUDA "OFF")
set(HAVE_TPETRA_INST_CUDA "OFF")
set(Tpetra_INST_HIP "OFF")
set(HAVE_TPETRA_INST_HIP "OFF")
set(Tpetra_INST_SYCL "OFF")
set(HAVE_TPETRA_INST_SYCL "OFF")
set(Tpetra_ALLOCATE_IN_SHARED_SPACE "OFF")
set(HAVE_TPETRA_SHARED_ALLOCS "OFF")
set(Tpetra_ASSUME_CUDA_AWARE_MPI "OFF")
set(TPETRA_ASSUME_CUDA_AWARE_MPI "OFF")
set(Tpetra_ASSUME_GPU_AWARE_MPI "OFF")
set(TPETRA_ASSUME_GPU_AWARE_MPI "OFF")
set(Tpetra_INST_INT_INT "ON")
set(HAVE_TPETRA_INST_INT_INT "ON")
set(Tpetra_INST_INT_UNSIGNED "OFF")
set(HAVE_TPETRA_INST_INT_UNSIGNED "OFF")
set(Tpetra_INST_INT_LONG_LONG "OFF")
set(HAVE_TPETRA_INST_INT_LONG_LONG "OFF")
set(Tpetra_INST_INT_LONG "OFF")
set(HAVE_TPETRA_INST_INT_LONG "OFF")
set(Tpetra_INST_INT_UNSIGNED_LONG "OFF")
set(HAVE_TPETRA_INST_INT_UNSIGNED_LONG "OFF")
set(Tpetra_INST_DOUBLE "ON")
set(HAVE_TPETRA_INST_DOUBLE "ON")
set(Tpetra_INST_FLOAT "OFF")
set(HAVE_TPETRA_INST_FLOAT "OFF")
set(Tpetra_INST_LONG_DOUBLE "OFF")
set(HAVE_TPETRA_INST_LONG_DOUBLE "OFF")
set(Tpetra_INST_COMPLEX_DOUBLE "ON")
set(HAVE_TPETRA_INST_COMPLEX_DOUBLE "ON")
set(Tpetra_INST_COMPLEX_FLOAT "OFF")
set(HAVE_TPETRA_INST_COMPLEX_FLOAT "OFF")
set(Tpetra_INST_FLOAT128 "OFF")
set(HAVE_TPETRA_INST_FLOAT128 "OFF")
set(Tpetra_INST_QD_REAL "OFF")
set(HAVE_TPETRA_INST_QD_REAL "OFF")
set(Tpetra_INST_DD_REAL "OFF")
set(HAVE_TPETRA_INST_DD_REAL "OFF")
set(Tpetra_ENABLE_Experimental "NO")
set(HAVE_TPETRA_EXPERIMENTAL "OFF")

# Include configuration of dependent packages
if (NOT TARGET TpetraCore::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TpetraCore/TpetraCoreConfig.cmake")
endif()
if (NOT TARGET TpetraTSQR::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TpetraTSQR/TpetraTSQRConfig.cmake")
endif()

# Import Tpetra targets
include("${CMAKE_CURRENT_LIST_DIR}/TpetraTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Tpetra_IS_TRIBITS_COMPLIANT TRUE)
set(Tpetra_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Tpetra_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Tpetra_EXPORTED_PACKAGE_LIBS_NAMES "kokkostsqr;tpetraclassic;tpetra;tpetrainout;tpetraext")

foreach(libname IN LISTS Tpetra_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Tpetra::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Tpetra::${libname}', or better yet,"
      " 'Tpetra::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Tpetra'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Tpetra_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
