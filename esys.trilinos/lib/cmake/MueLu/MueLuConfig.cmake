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
# CMake variable for use by Trilinos/MueLu clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "MueLu requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/MueLu build
## ---------------------------------------------------------------------------

set(MueLu_CXX_COMPILER "/usr/bin/g++")

set(MueLu_C_COMPILER "/usr/bin/gcc")

set(MueLu_Fortran_COMPILER "")
# Deprecated!
set(MueLu_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/MueLu build
## ---------------------------------------------------------------------------

## Give the build type
set(MueLu_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(MueLu_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(MueLu_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(MueLu_Fortran_FLAGS [[ ]])
# Deprecated
set(MueLu_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(MueLu_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(MueLu_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(MueLu_BUILD_SHARED_LIBS "ON")

set(MueLu_LINKER /usr/bin/ld)
set(MueLu_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(MueLu_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(MueLu_LIBRARIES MueLu::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(MueLu_MPI_LIBRARIES "")
set(MueLu_MPI_LIBRARY_DIRS "")
set(MueLu_MPI_INCLUDE_DIRS "")
set(MueLu_MPI_EXEC "")
set(MueLu_MPI_EXEC_MAX_NUMPROCS "")
set(MueLu_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(MueLu_ENABLE_Teuchos ON)
set(MueLu_ENABLE_Tpetra ON)
set(MueLu_ENABLE_Xpetra ON)
set(MueLu_ENABLE_Kokkos ON)
set(MueLu_ENABLE_KokkosKernels ON)
set(MueLu_ENABLE_BLAS ON)
set(MueLu_ENABLE_LAPACK ON)
set(MueLu_ENABLE_Amesos ON)
set(MueLu_ENABLE_Amesos2 ON)
set(MueLu_ENABLE_Belos ON)
set(MueLu_ENABLE_Epetra ON)
set(MueLu_ENABLE_EpetraExt ON)
set(MueLu_ENABLE_Teko OFF)
set(MueLu_ENABLE_Ifpack ON)
set(MueLu_ENABLE_Ifpack2 ON)
set(MueLu_ENABLE_Intrepid2 OFF)
set(MueLu_ENABLE_ML ON)
set(MueLu_ENABLE_Zoltan OFF)
set(MueLu_ENABLE_Zoltan2Core OFF)
set(MueLu_ENABLE_Stratimikos OFF)
set(MueLu_ENABLE_Thyra OFF)
set(MueLu_ENABLE_ThyraTpetraAdapters OFF)
set(MueLu_ENABLE_Isorropia OFF)
set(MueLu_ENABLE_Boost ON)
set(MueLu_ENABLE_MATLAB OFF)
set(MueLu_ENABLE_AmgX OFF)
set(MueLu_ENABLE_ViennaCL OFF)
set(MueLu_ENABLE_MKL OFF)
set(MueLu_ENABLE_Avatar OFF)
set(MueLu_ENABLE_CUSPARSE OFF)
set(MueLu_ENABLE_MAGMASparse OFF)
set(MueLu_ENABLE_mlpack OFF)

# Exported cache variables
set(MueLu_ENABLE_EXPLICIT_INSTANTIATION "ON")
set(HAVE_MUELU_EXPLICIT_INSTANTIATION "ON")
set(MueLu_ENABLE_DEBUG "OFF")
set(HAVE_MUELU_DEBUG "OFF")
set(MueLu_ENABLE_DEPRECATED_CODE "YES")
set(HAVE_MUELU_DEPRECATED_CODE "ON")
set(MueLu_ENABLE_DEPRECATED_TESTS "NO")
set(HAVE_MUELU_DEPRECATED_TESTS "OFF")
set(MueLu_ENABLE_Experimental "NO")
set(HAVE_MUELU_EXPERIMENTAL "OFF")
set(MueLu_ENABLE_ADDITIVE_VARIANT "NO")
set(HAVE_MUELU_ADDITIVE_VARIANT "OFF")
set(MueLu_ENABLE_REGION_SPLITTING "NO")
set(HAVE_MUELU_REGION_SPLITTING "OFF")
set(MueLu_ENABLE_Tutorial "MueLu_ENABLE_Tutorial_DEFAULT")
set(MueLu_ENABLE_SPLIT_ETI_CPP_FILES "NO")
set(HAVE_MUELU_SPLIT_ETI_CPP_FILES "OFF")
set(MueLu_ENABLE_GOOGLE_PERFTOOLS "OFF")
set(HAVE_MUELU_GOOGLE_PERFTOOLS "OFF")
set(MueLu_ENABLE_Boost_for_real "OFF")
set(HAVE_MUELU_BOOST_FOR_REAL "OFF")
set(MueLu_ENABLE_ML_MMM "OFF")
set(HAVE_MUELU_ML_MMM "OFF")

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Tpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Tpetra/TpetraConfig.cmake")
endif()
if (NOT TARGET Xpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Xpetra/XpetraConfig.cmake")
endif()
if (NOT TARGET Kokkos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Kokkos/KokkosConfig.cmake")
endif()
if (NOT TARGET KokkosKernels::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../KokkosKernels/KokkosKernelsConfig.cmake")
endif()
if (NOT TARGET BLAS::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/BLAS/BLASConfig.cmake")
endif()
if (NOT TARGET LAPACK::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/LAPACK/LAPACKConfig.cmake")
endif()
if (NOT TARGET Amesos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Amesos/AmesosConfig.cmake")
endif()
if (NOT TARGET Amesos2::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Amesos2/Amesos2Config.cmake")
endif()
if (NOT TARGET Belos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Belos/BelosConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET EpetraExt::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../EpetraExt/EpetraExtConfig.cmake")
endif()
if (NOT TARGET Ifpack::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Ifpack/IfpackConfig.cmake")
endif()
if (NOT TARGET Ifpack2::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Ifpack2/Ifpack2Config.cmake")
endif()
if (NOT TARGET ML::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../ML/MLConfig.cmake")
endif()
if (NOT TARGET Boost::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../../external_packages/Boost/BoostConfig.cmake")
endif()

# Import MueLu targets
include("${CMAKE_CURRENT_LIST_DIR}/MueLuTargets.cmake")

# Standard TriBITS-compliant external package variables
set(MueLu_IS_TRIBITS_COMPLIANT TRUE)
set(MueLu_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(MueLu_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(MueLu_EXPORTED_PACKAGE_LIBS_NAMES "muelu;muelu-adapters")

foreach(libname IN LISTS MueLu_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE MueLu::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'MueLu::${libname}', or better yet,"
      " 'MueLu::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'MueLu'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'MueLu_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
