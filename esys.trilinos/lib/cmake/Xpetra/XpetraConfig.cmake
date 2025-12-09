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
# CMake variable for use by Trilinos/Xpetra clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Xpetra requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Xpetra build
## ---------------------------------------------------------------------------

set(Xpetra_CXX_COMPILER "/usr/bin/g++")

set(Xpetra_C_COMPILER "/usr/bin/gcc")

set(Xpetra_Fortran_COMPILER "")
# Deprecated!
set(Xpetra_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Xpetra build
## ---------------------------------------------------------------------------

## Give the build type
set(Xpetra_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Xpetra_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Xpetra_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(Xpetra_Fortran_FLAGS [[ ]])
# Deprecated
set(Xpetra_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(Xpetra_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Xpetra_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(Xpetra_BUILD_SHARED_LIBS "ON")

set(Xpetra_LINKER /usr/bin/ld)
set(Xpetra_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Xpetra_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(Xpetra_LIBRARIES Xpetra::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Xpetra_MPI_LIBRARIES "")
set(Xpetra_MPI_LIBRARY_DIRS "")
set(Xpetra_MPI_INCLUDE_DIRS "")
set(Xpetra_MPI_EXEC "")
set(Xpetra_MPI_EXEC_MAX_NUMPROCS "")
set(Xpetra_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Xpetra_ENABLE_Teuchos ON)
set(Xpetra_ENABLE_Kokkos ON)
set(Xpetra_ENABLE_BLAS ON)
set(Xpetra_ENABLE_LAPACK ON)
set(Xpetra_ENABLE_Epetra ON)
set(Xpetra_ENABLE_EpetraExt ON)
set(Xpetra_ENABLE_Tpetra ON)
set(Xpetra_ENABLE_KokkosKernels ON)
set(Xpetra_ENABLE_Thyra OFF)

# Exported cache variables
set(Xpetra_ENABLE_DEBUG "OFF")
set(HAVE_XPETRA_DEBUG "OFF")
set(Xpetra_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_XPETRA_EXPLICIT_INSTANTIATION "ON")
set(Xpetra_ENABLE_DEPRECATED_CODE "ON")
set(XPETRA_ENABLE_DEPRECATED_CODE "ON")
set(Xpetra_ENABLE_Experimental "NO")
set(HAVE_XPETRA_EXPERIMENTAL "OFF")
set(Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES "OFF")
set(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES "OFF")
set(Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES "ON")
set(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES "ON")
set(XPETRA_ENABLE_PROFILING "OFF")
set(HAVE_XPETRA_PROFILING "OFF")

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Kokkos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Kokkos/KokkosConfig.cmake")
endif()
if (NOT TARGET BLAS::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/BLAS/BLASConfig.cmake")
endif()
if (NOT TARGET LAPACK::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/LAPACK/LAPACKConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET EpetraExt::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../EpetraExt/EpetraExtConfig.cmake")
endif()
if (NOT TARGET Tpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Tpetra/TpetraConfig.cmake")
endif()
if (NOT TARGET KokkosKernels::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../KokkosKernels/KokkosKernelsConfig.cmake")
endif()

# Import Xpetra targets
include("${CMAKE_CURRENT_LIST_DIR}/XpetraTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Xpetra_IS_TRIBITS_COMPLIANT TRUE)
set(Xpetra_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Xpetra_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Xpetra_EXPORTED_PACKAGE_LIBS_NAMES "xpetra;xpetra-sup")

foreach(libname IN LISTS Xpetra_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Xpetra::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Xpetra::${libname}', or better yet,"
      " 'Xpetra::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Xpetra'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Xpetra_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
