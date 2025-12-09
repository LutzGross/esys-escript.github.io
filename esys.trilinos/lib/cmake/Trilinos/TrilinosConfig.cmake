# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
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
# CMake variable for use by Trilinos clients. 
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

if(CMAKE_VERSION VERSION_LESS 3.3)
  set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    "Trilinos requires CMake 3.3 or later for 'if (... IN_LIST ...)'"
    )
  set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  return()
endif()
cmake_minimum_required(VERSION 3.3...3.23.0)

## ---------------------------------------------------------------------------
## Compilers used by Trilinos build
## ---------------------------------------------------------------------------

set(Trilinos_CXX_COMPILER "/usr/bin/g++")

set(Trilinos_C_COMPILER "/usr/bin/gcc")

set(Trilinos_Fortran_COMPILER "")

## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos build
## ---------------------------------------------------------------------------

set(Trilinos_CMAKE_BUILD_TYPE "RELEASE")

set(Trilinos_CXX_COMPILER_FLAGS [[   -Wno-unused-parameter -Wno-error=implicit-function-declaration   -O3 -DNDEBUG]])

set(Trilinos_C_COMPILER_FLAGS [[   -Wno-error=implicit-function-declaration -fopenmp -O3 -DNDEBUG]])

set(Trilinos_Fortran_COMPILER_FLAGS [[ -fopenmp ]])

## Extra link flags (e.g., specification of fortran libraries)
set(Trilinos_EXTRA_LD_FLAGS [[]])

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty. 
set(Trilinos_SHARED_LIB_RPATH_COMMAND "-Wl,-rpath,/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/lib")
set(Trilinos_BUILD_SHARED_LIBS "ON")

set(Trilinos_LINKER /usr/bin/ld)
set(Trilinos_AR /usr/bin/ar)


## ---------------------------------------------------------------------------
## Set library specifications and paths 
## ---------------------------------------------------------------------------

## The project version number
set(Trilinos_VERSION "16.0.0")

# For best practices in handling of components, see
# <http://www.cmake.org/cmake/help/v3.2/manual/cmake-developer.7.html#find-modules>.
#
# If components were requested, include only those. If not, include all of
# Trilinos.
if (Trilinos_FIND_COMPONENTS)
  set(COMPONENTS_LIST ${Trilinos_FIND_COMPONENTS})
else()
  set(COMPONENTS_LIST MueLu;Ifpack2;Amesos2;Belos;ML;Ifpack;Amesos;Galeri;AztecOO;Xpetra;TrilinosSS;Tpetra;TpetraCore;TpetraTSQR;EpetraExt;Triutils;Epetra;KokkosKernels;Teuchos;TeuchosKokkosComm;TeuchosKokkosCompat;TeuchosRemainder;TeuchosNumerics;TeuchosComm;TeuchosParameterList;TeuchosParser;TeuchosCore;Kokkos)
endif()

# Initialize Trilinos_FOUND with true, and set it to FALSE if any of
# the required components wasn't found.
set(Trilinos_FOUND TRUE)
set(Trilinos_NOT_FOUND_MESSAGE "")
set(selectedComponentsFound "")
foreach (comp IN ITEMS ${COMPONENTS_LIST})
  set(compPkgConfigFile
    ${CMAKE_CURRENT_LIST_DIR}/../${comp}/${comp}Config.cmake
    )
  if (EXISTS ${compPkgConfigFile})
    # Set Trilinos_<component>_FOUND.
    set(Trilinos_${comp}_FOUND TRUE)
    # Include the package file.
    include(${compPkgConfigFile})
    # Add variables to lists.
    list(APPEND Trilinos_LIBRARIES ${${comp}_LIBRARIES})
    list(APPEND selectedComponentsFound ${comp})
  else()
    set(Trilinos_${comp}_FOUND FALSE)
    if(Trilinos_FIND_REQUIRED_${comp})
      string(APPEND Trilinos_NOT_FOUND_MESSAGE
        "ERROR: Could not find component '${comp}'!\n")
      set(Trilinos_FOUND FALSE)
    endif()
  endif()
endforeach()

# Deprecated (see #299)!
set(Trilinos_INCLUDE_DIRS "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos/include")

# Remove duplicates in Trilinos_LIBRARIES
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

set(Trilinos_INSTALL_DIR "/home/lgross/PycharmProjects/esys-escript.github.io/esys.trilinos")
set(Trilinos_MPI_LIBRARIES "")
set(Trilinos_MPI_LIBRARY_DIRS "")
set(Trilinos_MPI_INCLUDE_DIRS "")
set(Trilinos_MPI_EXEC "")
set(Trilinos_MPI_EXEC_PRE_NUMPROCS_FLAGS "")
set(Trilinos_MPI_EXEC_MAX_NUMPROCS "")
set(Trilinos_MPI_EXEC_POST_NUMPROCS_FLAGS "")
set(Trilinos_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Compiler vendor identifications
## ---------------------------------------------------------------------------
set(Trilinos_SYSTEM_NAME "Linux")
set(Trilinos_CXX_COMPILER_ID "GNU")
set(Trilinos_C_COMPILER_ID "GNU")
set(Trilinos_Fortran_COMPILER_ID "")
set(Trilinos_Fortran_IMPLICIT_LINK_LIBRARIES "")

## ---------------------------------------------------------------------------
## Set useful general variables 
## ---------------------------------------------------------------------------

## The packages enabled for this project
set(Trilinos_PACKAGE_LIST "MueLu;Ifpack2;Amesos2;Belos;ML;Ifpack;Amesos;Galeri;AztecOO;Xpetra;TrilinosSS;Tpetra;TpetraCore;TpetraTSQR;EpetraExt;Triutils;Epetra;KokkosKernels;Teuchos;TeuchosKokkosComm;TeuchosKokkosCompat;TeuchosRemainder;TeuchosNumerics;TeuchosComm;TeuchosParameterList;TeuchosParser;TeuchosCore;Kokkos")

## The selected packages for this project
set(Trilinos_SELECTED_PACKAGE_LIST "${selectedComponentsFound}")

## ---------------------------------------------------------------------------
## Modern CMake (IMPORTED) targets
## ---------------------------------------------------------------------------

# Trilinos::all_libs  (Does *not* depend on COMPONENTS)
if (NOT TARGET Trilinos::all_libs)
  set(Trilinos_ALL_PACKAGES_TARGETS)
  foreach (pkg IN ITEMS MueLu;Ifpack2;Amesos2;Belos;ML;Ifpack;Amesos;Galeri;AztecOO;Xpetra;TrilinosSS;Tpetra;TpetraCore;TpetraTSQR;EpetraExt;Triutils;Epetra;KokkosKernels;Teuchos;TeuchosKokkosComm;TeuchosKokkosCompat;TeuchosRemainder;TeuchosNumerics;TeuchosComm;TeuchosParameterList;TeuchosParser;TeuchosCore;Kokkos)
    list(APPEND Trilinos_ALL_PACKAGES_TARGETS ${pkg}::all_libs)
  endforeach()
  add_library(Trilinos::all_libs IMPORTED INTERFACE GLOBAL)
  target_link_libraries(Trilinos::all_libs
  INTERFACE ${Trilinos_ALL_PACKAGES_TARGETS} )
endif()

# Trilinos::all_selected_libs  (Depend on COMPONENTS)
if (NOT TARGET Trilinos::all_selected_libs)
  set(Trilinos_ALL_SELECTED_PACKAGES_TARGETS)
  foreach (pkg IN ITEMS ${selectedComponentsFound})
    list(APPEND Trilinos_ALL_SELECTED_PACKAGES_TARGETS ${pkg}::all_libs)
  endforeach()
  add_library(Trilinos::all_selected_libs IMPORTED INTERFACE GLOBAL)
  target_link_libraries(Trilinos::all_selected_libs
    INTERFACE ${Trilinos_ALL_SELECTED_PACKAGES_TARGETS} )
endif()
