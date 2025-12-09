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
# CMake variable for use by Trilinos/TeuchosCore clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "TeuchosCore requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/TeuchosCore build
## ---------------------------------------------------------------------------

set(TeuchosCore_CXX_COMPILER "/usr/bin/g++")

set(TeuchosCore_C_COMPILER "/usr/bin/gcc")

set(TeuchosCore_Fortran_COMPILER "")
# Deprecated!
set(TeuchosCore_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/TeuchosCore build
## ---------------------------------------------------------------------------

## Give the build type
set(TeuchosCore_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(TeuchosCore_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(TeuchosCore_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(TeuchosCore_Fortran_FLAGS [[ ]])
# Deprecated
set(TeuchosCore_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(TeuchosCore_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(TeuchosCore_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(TeuchosCore_BUILD_SHARED_LIBS "ON")

set(TeuchosCore_LINKER /usr/bin/ld)
set(TeuchosCore_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(TeuchosCore_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(TeuchosCore_LIBRARIES TeuchosCore::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(TeuchosCore_MPI_LIBRARIES "")
set(TeuchosCore_MPI_LIBRARY_DIRS "")
set(TeuchosCore_MPI_INCLUDE_DIRS "")
set(TeuchosCore_MPI_EXEC "")
set(TeuchosCore_MPI_EXEC_MAX_NUMPROCS "")
set(TeuchosCore_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(TeuchosCore_ENABLE_Kokkos ON)
set(TeuchosCore_ENABLE_BinUtils OFF)
set(TeuchosCore_ENABLE_Boost ON)
set(TeuchosCore_ENABLE_MPI OFF)
set(TeuchosCore_ENABLE_ARPREC OFF)
set(TeuchosCore_ENABLE_QD OFF)
set(TeuchosCore_ENABLE_QT OFF)
set(TeuchosCore_ENABLE_quadmath OFF)
set(TeuchosCore_ENABLE_Pthread OFF)
set(TeuchosCore_ENABLE_Valgrind OFF)

# Exported cache variables
set(Teuchos_ENABLE_DEBUG "OFF")
set(HAVE_TEUCHOS_DEBUG "OFF")
set(Teuchos_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_TEUCHOS_EXPLICIT_INSTANTIATION "ON")
set(Teuchos_ENABLE_LONG_DOUBLE "FALSE")
set(HAVE_TEUCHOS_LONG_DOUBLE "OFF")
set(Teuchos_ENABLE_FLOAT "FALSE")
set(HAVE_TEUCHOS_FLOAT "OFF")
set(Teuchos_ENABLE_COMPLEX "ON")
set(HAVE_TEUCHOS_COMPLEX "ON")
set(Teuchos_INST_FLOAT "FALSE")
set(HAVE_TEUCHOS_INST_FLOAT "OFF")
set(Teuchos_INST_COMPLEX_FLOAT "OFF")
set(HAVE_TEUCHOS_INST_COMPLEX_FLOAT "OFF")
set(Teuchos_INST_COMPLEX_DOUBLE "TRUE")
set(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE "ON")
set(Teuchos_ENABLE_THREAD_SAFE "ON")
set(HAVE_TEUCHOS_THREAD_SAFE "ON")
set(Teuchos_ENABLE_ABC "OFF")
set(HAVE_TEUCHOS_ARRAY_BOUNDSCHECK "OFF")
set(Teuchos_MODIFY_DEFAULTS_DURING_VALIDATION "ON")
set(HAVE_TEUCHOS_MODIFY_DEFAULTS_DURING_VALIDATION "ON")
set(Teuchos_ENABLE_DEBUG_RCP_NODE_TRACING "OFF")
set(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING "OFF")
set(Teuchos_ENABLE_EXTENDED "ON")
set(HAVE_TEUCHOS_EXTENDED "ON")
set(Teuchos_ENABLE_C_EXCEPTIONS "OFF")
set(HAVE_TEUCHOS_C_EXCEPTIONS "OFF")
set(Teuchos_ENABLE_GCC_DEMANGLE "ON")
set(HAVE_TEUCHOS_DEMANGLE "ON")
set(Teuchos_ENABLE_STACKED_TIMER_IN_TIME_MONITOR "ON")
set(HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER "ON")
set(Teuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS "OFF")
set(HAVE_TEUCHOS_GLOBALLY_REDUCE_UNITTEST_RESULTS "OFF")
set(Teuchos_KOKKOS_PROFILING "ON")
set(HAVE_TEUCHOS_KOKKOS_PROFILING "ON")
set(Teuchos_TIMER_KOKKOS_FENCE "OFF")
set(HAVE_TEUCHOS_TIMER_KOKKOS_FENCE "OFF")

# Include configuration of dependent packages
if (NOT TARGET Kokkos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Kokkos/KokkosConfig.cmake")
endif()
if (NOT TARGET Boost::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/Boost/BoostConfig.cmake")
endif()

# Import TeuchosCore targets
include("${CMAKE_CURRENT_LIST_DIR}/TeuchosCoreTargets.cmake")

# Standard TriBITS-compliant external package variables
set(TeuchosCore_IS_TRIBITS_COMPLIANT TRUE)
set(TeuchosCore_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(TeuchosCore_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(TeuchosCore_EXPORTED_PACKAGE_LIBS_NAMES "teuchoscore")

foreach(libname IN LISTS TeuchosCore_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE TeuchosCore::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'TeuchosCore::${libname}', or better yet,"
      " 'TeuchosCore::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'TeuchosCore'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'TeuchosCore_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
