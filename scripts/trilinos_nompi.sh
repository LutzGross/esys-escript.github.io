#!/bin/bash

rm -f CMakeCache.txt

TRI_INSTALL_PREFIX=$1
# $2 = C compiler
# $3 = C++ compiler
# $4 = OpenMP flag (not used, controlled by Trilinos_ENABLE_OpenMP)
# $5 = Trilinos source directory
# $6 = METIS enable (ON/OFF)
# $7 = METIS include path
# $8 = METIS library path
# $9 = ParMETIS enable (ignored for nompi build)
# ${10} = ParMETIS include path (ignored)
# ${11} = ParMETIS library path (ignored)

METIS_ENABLE=$6
METIS_INC=$7
METIS_LIB=$8

# Build METIS cmake options
METIS_OPTS="-D TPL_ENABLE_METIS=${METIS_ENABLE}"
if [ "${METIS_ENABLE}" = "ON" ] && [ -n "${METIS_INC}" ]; then
    METIS_OPTS="${METIS_OPTS} -D TPL_METIS_INCLUDE_DIRS=${METIS_INC}"
    METIS_OPTS="${METIS_OPTS} -D TPL_METIS_LIBRARIES=${METIS_LIB}/libmetis.so"
fi

cmake \
      -D CMAKE_INSTALL_PREFIX=$TRI_INSTALL_PREFIX\
      -D CMAKE_CXX_STANDARD=20 \
      -D Trilinos_ENABLE_Fortran=OFF \
      -D CMAKE_C_COMPILER=$2 \
      -D CMAKE_C_FLAGS=" -Wno-error=implicit-function-declaration " \
      -D CMAKE_CXX_COMPILER=$3 \
      -D CMAKE_CXX_FLAGS=" -Wno-unused-parameter -Wno-error=implicit-function-declaration " \
      -D BUILD_SHARED_LIBS=ON \
      -D TPL_ENABLE_BLAS=ON \
      -D TPL_ENABLE_Boost=ON \
      -D TPL_ENABLE_Cholmod=OFF \
      -D TPL_ENABLE_SCALAPACK=OFF \
      ${METIS_OPTS} \
      -D TPL_ENABLE_ParMETIS=OFF \
      -D TPL_ENABLE_MUMPS=OFF \
      -D TPL_ENABLE_SuperLU=OFF \
      -D TPL_ENABLE_UMFPACK=OFF \
      -D Trilinos_ENABLE_Amesos2=ON \
      -D Trilinos_ENABLE_Belos=ON \
      -D Trilinos_ENABLE_Ifpack2=ON \
      -D Trilinos_ENABLE_Kokkos=ON \
      -D Trilinos_ENABLE_MueLu=ON \
      -D Trilinos_ENABLE_Teuchos=ON \
      -D Trilinos_ENABLE_Tpetra=ON \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
      -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON \
      -D Tpetra_INST_COMPLEX_DOUBLE=ON \
      -D Trilinos_ENABLE_COMPLEX_DOUBLE=ON \
      -D Teuchos_ENABLE_COMPLEX=ON \
      -D Teuchos_ENABLE_THREAD_SAFE=ON \
      -D Tpetra_INST_INT_INT=ON \
      -D Trilinos_ENABLE_OpenMP=ON \
      -D TPL_ENABLE_MPI=OFF \
      -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -D Kokkos_ENABLE_COMPILER_WARNINGS=ON \
      -D Amesos2_ENABLE_Basker=ON \
      -D Tpetra_INST_SERIAL:BOOL=ON \
      -D Trilinos_ENABLE_TESTS=OFF \
      -D Trilinos_ENABLE_EXAMPLES=OFF \
      $5
