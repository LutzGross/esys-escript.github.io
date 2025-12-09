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
# CMake variable for use by Trilinos/Galeri clients.
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Galeri requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Galeri build
## ---------------------------------------------------------------------------

set(Galeri_CXX_COMPILER "/usr/bin/g++")

set(Galeri_C_COMPILER "/usr/bin/gcc")

set(Galeri_Fortran_COMPILER "")
# Deprecated!
set(Galeri_FORTRAN_COMPILER "") 


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Galeri build
## ---------------------------------------------------------------------------

## Give the build type
set(Galeri_CMAKE_BUILD_TYPE "RELEASE")

## Set compiler flags, including those determined by build type
set(Galeri_CXX_FLAGS [[ -O3 -DNDEBUG]])

set(Galeri_C_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(Galeri_Fortran_FLAGS [[ ]])
# Deprecated
set(Galeri_FORTRAN_FLAGS [[ ]])

## Extra link flags (e.g., specification of fortran libraries)
set(Galeri_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty.
set(Galeri_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(Galeri_BUILD_SHARED_LIBS "ON")

set(Galeri_LINKER /usr/bin/ld)
set(Galeri_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths
## ---------------------------------------------------------------------------

## Base install location (if not in the build tree)
set(Galeri_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")

## List of package libraries
set(Galeri_LIBRARIES Galeri::all_libs)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Galeri_MPI_LIBRARIES "")
set(Galeri_MPI_LIBRARY_DIRS "")
set(Galeri_MPI_INCLUDE_DIRS "")
set(Galeri_MPI_EXEC "")
set(Galeri_MPI_EXEC_MAX_NUMPROCS "")
set(Galeri_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables
## ---------------------------------------------------------------------------

# Enables/Disables for upstream package dependencies
set(Galeri_ENABLE_Teuchos ON)
set(Galeri_ENABLE_Epetra ON)
set(Galeri_ENABLE_EpetraExt ON)
set(Galeri_ENABLE_Xpetra ON)
set(Galeri_ENABLE_Tpetra ON)
set(Galeri_ENABLE_Triutils ON)

# Exported cache variables

# Include configuration of dependent packages
if (NOT TARGET Teuchos::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Teuchos/TeuchosConfig.cmake")
endif()
if (NOT TARGET Epetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Epetra/EpetraConfig.cmake")
endif()
if (NOT TARGET EpetraExt::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../EpetraExt/EpetraExtConfig.cmake")
endif()
if (NOT TARGET Xpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Xpetra/XpetraConfig.cmake")
endif()
if (NOT TARGET Tpetra::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Tpetra/TpetraConfig.cmake")
endif()
if (NOT TARGET Triutils::all_libs)
  include("${CMAKE_CURRENT_LIST_DIR}/../Triutils/TriutilsConfig.cmake")
endif()

# Import Galeri targets
include("${CMAKE_CURRENT_LIST_DIR}/GaleriTargets.cmake")

# Standard TriBITS-compliant external package variables
set(Galeri_IS_TRIBITS_COMPLIANT TRUE)
set(Galeri_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")
set(Galeri_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR "${CMAKE_CURRENT_LIST_DIR}")


## ----------------------------------------------------------------------------
## Create deprecated non-namespaced library targets for backwards compatibility
## ----------------------------------------------------------------------------

set(Galeri_EXPORTED_PACKAGE_LIBS_NAMES "galeri-epetra;galeri-xpetra")

foreach(libname IN LISTS Galeri_EXPORTED_PACKAGE_LIBS_NAMES)
  if (NOT TARGET ${libname})
    add_library(${libname} INTERFACE IMPORTED)
    target_link_libraries(${libname}
       INTERFACE Galeri::${libname})
    set(deprecationMessage
      "WARNING: The non-namespaced target '${libname}' is deprecated!"
      "  If always using newer versions of the project 'Trilinos', then use the"
      " new namespaced target 'Galeri::${libname}', or better yet,"
      " 'Galeri::all_libs' to be less sensitive to changes in the definition"
      " of targets in the package 'Galeri'.  Or, to maintain compatibility with"
      " older or newer versions the project 'Trilinos', instead link against the"
      " libraries specified by the variable 'Galeri_LIBRARIES'."
      )
    string(REPLACE ";" "" deprecationMessage "${deprecationMessage}")
    set_target_properties(${libname}
      PROPERTIES DEPRECATION "${deprecationMessage}" )
  endif()
endforeach()
