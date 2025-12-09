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
# CMake variable for use by Trilinos/TpetraCore clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "TpetraCore requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/TpetraCore build
## ---------------------------------------------------------------------------

set(TpetraCore_CXX_COMPILER "/usr/bin/g++")

set(TpetraCore_C_COMPILER "/usr/bin/gcc")

set(TpetraCore_Fortran_COMPILER "")
# Deprecated!
set(TpetraCore_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/TpetraCore build
## ---------------------------------------------------------------------------

## Give the build type
set(TpetraCore_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(TpetraCore_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(TpetraCore_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(TpetraCore_Fortran_FLAGS [[ ]])
# Deprecated
set(TpetraCore_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(TpetraCore_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(TpetraCore_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(TpetraCore_BUILD_SHARED_LIBS "ON")

set(TpetraCore_LINKER /usr/bin/ld)
set(TpetraCore_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(TpetraCore_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(TpetraCore_LIBRARIES TpetraCore::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(TpetraCore_MPI_LIBRARIES "")
set(TpetraCore_MPI_LIBRARY_DIRS "")
set(TpetraCore_MPI_INCLUDE_DIRS "")
set(TpetraCore_MPI_EXEC "")
set(TpetraCore_MPI_EXEC_MAX_NUMPROCS "")
set(TpetraCore_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(TpetraCore_ENABLE_Teuchos ON)
set(TpetraCore_ENABLE_Kokkos ON)
set(TpetraCore_ENABLE_TeuchosKokkosCompat ON)
set(TpetraCore_ENABLE_TeuchosKokkosComm ON)
set(TpetraCore_ENABLE_KokkosKernels ON)
set(TpetraCore_ENABLE_Epetra ON)
set(TpetraCore_ENABLE_TpetraTSQR ON)
set(TpetraCore_ENABLE_TeuchosNumerics ON)
set(TpetraCore_ENABLE_MPI OFF)
set(TpetraCore_ENABLE_CUDA OFF)
set(TpetraCore_ENABLE_QD OFF)
set(TpetraCore_ENABLE_quadmath OFF)
set(TpetraCore_ENABLE_mpi_advance OFF)

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
set(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_TPETRACORE_EXPLICIT_INSTANTIATION "ON")
set(Tpetra_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_TPETRA_EXPLICIT_INSTANTIATION "ON")
set(TpetraCore_ENABLE_SS_TESTING "OFF")
set(HAVE_TPETRA_ENABLE_SS_TESTING "OFF")
set(TpetraCore_ENABLE_TSQR "ON")
set(HAVE_TPETRA_TSQR "ON")
set(TpetraCore_ENABLE_BUGTESTS "ON")
set(HAVE_TPETRA_BUGTESTS "ON")
set(Tpetra_ENABLE_Reduced_ETI "OFF")
set(HAVE_TPETRA_REDUCED_ETI "OFF")
set(TpetraCore_Threaded_MKL "OFF")
set(HAVE_TPETRA_THREADED_MKL "OFF")
set(Tpetra_THROW_Warnings "OFF")
set(HAVE_TPETRA_THROW_WARNINGS "OFF")
set(Tpetra_PRINT_Warnings "OFF")
set(HAVE_TPETRA_PRINT_WARNINGS "OFF")
set(Tpetra_PRINT_Efficiency_Warnings "OFF")
set(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS "OFF")
set(Tpetra_THROW_Abuse_Warnings "OFF")
set(HAVE_TPETRA_THROW_ABUSE_WARNINGS "OFF")
set(Tpetra_PRINT_Abuse_Warnings "OFF")
set(HAVE_TPETRA_PRINT_ABUSE_WARNINGS "OFF")
set(Tpetra_USE_MURMUR_HASH "OFF")
set(TPETRA_USE_MURMUR_HASH "OFF")
set(Tpetra_ENABLE_MMM_Timings "OFF")
set(HAVE_TPETRA_MMM_TIMINGS "OFF")
set(Tpetra_ENABLE_Distributor_Timings "OFF")
set(HAVE_TPETRA_DISTRIBUTOR_TIMINGS "OFF")
set(Tpetra_ENABLE_MMM_Statistics "OFF")
set(HAVE_TPETRA_MMM_STATISTICS "OFF")
set(Tpetra_BCRS_Point_Import "OFF")
set(HAVE_TPETRA_BCRS_DO_POINT_IMPORT "OFF")

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Kokkos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Kokkos/KokkosConfig.cmake")
endif()
if (NOT TARGET TeuchosKokkosCompat::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosKokkosCompat/TeuchosKokkosCompatConfig.cmake")
endif()
if (NOT TARGET TeuchosKokkosComm::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosKokkosComm/TeuchosKokkosCommConfig.cmake")
endif()
if (NOT TARGET KokkosKernels::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../KokkosKernels/KokkosKernelsConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET TpetraTSQR::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TpetraTSQR/TpetraTSQRConfig.cmake")
endif()
if (NOT TARGET TeuchosNumerics::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TeuchosNumerics/TeuchosNumericsConfig.cmake")
endif()

# Import TpetraCore targets
include("${CMAKE_CURRENT_LIST_DIR}/TpetraCoreTargets.cmake")

# Standard TriBITS-compliant external package variables
set(TpetraCore_IS_TRIBITS_COMPLIANT TRUE)
set(TpetraCore_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(TpetraCore_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(TpetraCore_EXPORTED_PACKAGE_LIBS_NAMES "tpetraclassic;tpetra;tpetrainout;tpetraext")

foreach(libname IN LISTS TpetraCore_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE TpetraCore::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'TpetraCore::${libname}', or better yet,"
      " 'TpetraCore::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'TpetraCore'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'TpetraCore_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
