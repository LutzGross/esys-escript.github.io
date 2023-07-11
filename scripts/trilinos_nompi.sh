#!/bin/bash

rm -f CMakeCache.txt

TRI_INSTALL_PREFIX=$1"/escript_trilinos"

cmake \
      -D CMAKE_INSTALL_PREFIX=$TRI_INSTALL_PREFIX \
      -D Trilinos_ENABLE_CXX11=ON \
      -D Trilinos_ENABLE_Fortran=OFF \
      -D CMAKE_C_COMPILER=$2 \
      -D CMAKE_CXX_COMPILER=$3 \
      -D CMAKE_CXX_FLAGS=" -Wno-unused-parameter" \
      -D BUILD_SHARED_LIBS=ON \
      -D TPL_ENABLE_BLAS=ON \
      -D TPL_ENABLE_Boost=ON \
      -D TPL_ENABLE_Cholmod=OFF \
      -D TPL_ENABLE_SCALAPACK=OFF \
      -D SCALAPACK_LIBRARY_NAMES='libscalapack-openmpi.so' \
      -D TPL_ENABLE_ParMETIS=OFF \
      -D TPL_ENABLE_MUMPS=OFF \
      -D TPL_ENABLE_SuperLU=OFF \
      -D TPL_ENABLE_UMFPACK=OFF \
      -D Trilinos_ENABLE_Amesos=ON \
      -D Trilinos_ENABLE_Amesos2=ON \
      -D Trilinos_ENABLE_AztecOO=ON \
      -D Trilinos_ENABLE_Belos=ON \
      -D Trilinos_ENABLE_Ifpack=ON \
      -D Trilinos_ENABLE_Ifpack2=ON \
      -D Trilinos_ENABLE_Kokkos=ON \
      -D Trilinos_ENABLE_Komplex=ON \
      -D Trilinos_ENABLE_ML=ON \
      -D Trilinos_ENABLE_MueLu=ON \
      -D Trilinos_ENABLE_Teuchos=ON \
      -D Trilinos_ENABLE_Tpetra=ON \
      -D Trilinos_ENABLE_Epetra=ON \
      -D Trilinos_ENABLE_EpetraExt=ON \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
      -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON \
      -D Tpetra_INST_COMPLEX_DOUBLE=ON \
      -D Trilinos_ENABLE_COMPLEX_DOUBLE=ON \
      -D Teuchos_ENABLE_COMPLEX=ON \
      -D Tpetra_INST_INT_INT=ON \
      -D Tpetra_ENABLE_DEPRECATED_CODE=ON \
      -D Trilinos_ENABLE_OpenMP=ON \
      -D TPL_ENABLE_MPI=OFF \
      -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -D Kokkos_ENABLE_COMPILER_WARNINGS=ON \
      -D Amesos2_ENABLE_Basker=ON \
      -D Tpetra_INST_SERIAL:BOOL=ON \
      -D Trilinos_ENABLE_TESTS=OFF \
      ../trilinos_source
#-D TPL_Cholmod_LIBRARIES='libcholmod.so;libamd.so;libcolamd.so' \