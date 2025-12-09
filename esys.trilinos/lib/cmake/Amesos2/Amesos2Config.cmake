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
# CMake variable for use by Trilinos/Amesos2 clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Amesos2 requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Amesos2 build
## ---------------------------------------------------------------------------

set(Amesos2_CXX_COMPILER "/usr/bin/g++")

set(Amesos2_C_COMPILER "/usr/bin/gcc")

set(Amesos2_Fortran_COMPILER "")
# Deprecated!
set(Amesos2_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Amesos2 build
## ---------------------------------------------------------------------------

## Give the build type
set(Amesos2_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Amesos2_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Amesos2_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(Amesos2_Fortran_FLAGS [[ ]])
# Deprecated
set(Amesos2_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(Amesos2_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Amesos2_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(Amesos2_BUILD_SHARED_LIBS "ON")

set(Amesos2_LINKER /usr/bin/ld)
set(Amesos2_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Amesos2_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(Amesos2_LIBRARIES Amesos2::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Amesos2_MPI_LIBRARIES "")
set(Amesos2_MPI_LIBRARY_DIRS "")
set(Amesos2_MPI_INCLUDE_DIRS "")
set(Amesos2_MPI_EXEC "")
set(Amesos2_MPI_EXEC_MAX_NUMPROCS "")
set(Amesos2_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Amesos2_ENABLE_Teuchos ON)
set(Amesos2_ENABLE_Tpetra ON)
set(Amesos2_ENABLE_TrilinosSS ON)
set(Amesos2_ENABLE_Kokkos ON)
set(Amesos2_ENABLE_Epetra ON)
set(Amesos2_ENABLE_EpetraExt ON)
set(Amesos2_ENABLE_ShyLU_NodeBasker OFF)
set(Amesos2_ENABLE_ShyLU_NodeTacho OFF)
set(Amesos2_ENABLE_MPI OFF)
set(Amesos2_ENABLE_SuperLU OFF)
set(Amesos2_ENABLE_SuperLUMT OFF)
set(Amesos2_ENABLE_SuperLUDist OFF)
set(Amesos2_ENABLE_LAPACK ON)
set(Amesos2_ENABLE_UMFPACK OFF)
set(Amesos2_ENABLE_PARDISO_MKL OFF)
set(Amesos2_ENABLE_CSS_MKL OFF)
set(Amesos2_ENABLE_ParMETIS OFF)
set(Amesos2_ENABLE_METIS OFF)
set(Amesos2_ENABLE_Cholmod OFF)
set(Amesos2_ENABLE_MUMPS OFF)
set(Amesos2_ENABLE_STRUMPACK OFF)
set(Amesos2_ENABLE_CUSPARSE OFF)
set(Amesos2_ENABLE_CUSOLVER OFF)

# Exported cache variables
set(Amesos2_ENABLE_VERBOSE_DEBUG "OFF")
set(HAVE_AMESOS2_VERBOSE_DEBUG "OFF")
set(Amesos2_ENABLE_TIMERS "OFF")
set(HAVE_AMESOS2_TIMERS "OFF")
set(Amesos2_ENABLE_KLU2 "ON")
set(HAVE_AMESOS2_KLU2 "ON")
set(Amesos2_ENABLE_Basker "ON")
set(HAVE_AMESOS2_BASKER "ON")
set(Amesos2_ENABLE_Experimental "NO")
set(HAVE_AMESOS2_EXPERIMENTAL "OFF")
set(Amesos2_ENABLE_DEBUG "OFF")
set(HAVE_AMESOS2_DEBUG "OFF")
set(Amesos2_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_AMESOS2_EXPLICIT_INSTANTIATION "ON")

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Tpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Tpetra/TpetraConfig.cmake")
endif()
if (NOT TARGET TrilinosSS::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../TrilinosSS/TrilinosSSConfig.cmake")
endif()
if (NOT TARGET Kokkos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Kokkos/KokkosConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET EpetraExt::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../EpetraExt/EpetraExtConfig.cmake")
endif()
if (NOT TARGET LAPACK::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/LAPACK/LAPACKConfig.cmake")
endif()

# Import Amesos2 targets
include("${CMAKE_CURRENT_LIST_DIR}/Amesos2Targets.cmake")

# Standard TriBITS-compliant external package variables
set(Amesos2_IS_TRIBITS_COMPLIANT TRUE)
set(Amesos2_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Amesos2_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Amesos2_EXPORTED_PACKAGE_LIBS_NAMES "amesos2")

foreach(libname IN LISTS Amesos2_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Amesos2::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Amesos2::${libname}', or better yet,"
      " 'Amesos2::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Amesos2'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Amesos2_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
