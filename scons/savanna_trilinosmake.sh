#!/bin/bash
# this provides make file to build trilinos on savanna

rm -f CMakeCache.txt

TRI_INSTALL_PREFIX=$1"/escript_trilinos"

cmake \
      -D CMAKE_INSTALL_PREFIX=$TRI_INSTALL_PREFIX \
      -D Trilinos_ENABLE_CXX11=ON \
      -D CMAKE_CXX_COMPILER=mpiicpc \
      -D CMAKE_C_COMPILER=mpiicc \
      -D CMAKE_Fortran_COMPILER=mpiifort \
      -D BUILD_SHARED_LIBS=ON \
      -D TPL_ENABLE_BLAS=ON \
      -D TPL_ENABLE_Boost=ON \
      -D TPL_ENABLE_Cholmod=OFF \
      -D TPL_ENABLE_SCALAPACK=ON \
      -D SCALAPACK_LIBRARY_NAMES='libscalapack-openmpi.so' \
      -D TPL_ENABLE_ParMETIS=ON \
      -D TPL_ENABLE_MUMPS=ON \
      -D TPL_ENABLE_SuperLU=OFF \
      -D TPL_ENABLE_UMFPACK=OFF \
      -D Boost_INCLUDE_DIRS=/sw/libs/boost/1.61.0/include \
      -D TPL_BLACS_INCLUDE_DIRS=$MKLROOT/include \
      -D TPL_BLACS_LIBRARIES='mkl_scalapack_lp64;mkl_intel_lp64;mkl_core;mkl_intel_thread;mkl_blacs_intelmpi_lp64' \
      -D TPL_BLACS_LIBRARY_DIRS=$MKLROOT/lib/intel64 \
      -D TPL_BLAS_INCLUDE_DIRS=$MKLROOT/include \
      -D TPL_BLAS_LIBRARIES='mkl_scalapack_lp64;mkl_intel_lp64;mkl_core;mkl_intel_thread;mkl_blacs_intelmpi_lp64' \
      -D TPL_BLAS_LIBRARY_DIRS=$MKLROOT/lib/intel64 \
      -D TPL_LAPACK_INCLUDE_DIRS=$MKLROOT/include \
      -D TPL_LAPACK_LIBRARIES='mkl_scalapack_lp64;mkl_intel_lp64;mkl_core;mkl_intel_thread;mkl_blacs_intelmpi_lp64' \
      -D TPL_LAPACK_LIBRARY_DIRS=$MKLROOT/lib/intel64 \
      -D PARDISO_MKL_LIBRARY_DIRS=$MKLROOT/lib/intel64 \
      -D PARDISO_MKL_INCLUDE_DIRS=$MKLROOT/include \
      -D MKL_INCLUDE_DIRS=$MKLROOT/include \
      -D MKL_LIBRARY_DIRS=$MKLROOT/lib/intel64 \
      -D MUMPS_INCLUDE_DIRS=/sw/installer/MUMPS_5.1.2/include \
      -D MUMPS_LIBRARY_DIRS=/sw/installer/MUMPS_5.1.2/lib \
      -D MUMPS_LIBRARIES='cmumps;dmumps;smumps;zmumps;mumps_common;pord' \
      -D ParMETIS_INCLUDE_DIRS=/sw/libs/parmetis/4.0.3-impi/include \
      -D ParMETIS_LIBRARY_DIRS=/sw/libs/parmetis/4.0.3-impi/lib \
      -D TPL_SCALAPACK_INCLUDE_DIRS=$MKLROOT/include \
      -D TPL_SCALAPACK_LIBRARIES='mkl_scalapack_lp64;mkl_intel_lp64;mkl_core;mkl_intel_thread;mkl_blacs_intelmpi_lp64' \
      -D TPL_SCALAPACK_LIBRARY_DIRS=$MKLROOT/lib/intel64 \
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
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON \
      -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON \
      -D Tpetra_INST_COMPLEX_DOUBLE=ON \
      -D Trilinos_ENABLE_COMPLEX_DOUBLE=ON \
      -D Teuchos_ENABLE_COMPLEX=ON \
      -D Tpetra_INST_INT_INT=ON \
      -D Tpetra_ENABLE_DEPRECATED_CODE=ON \
      -D Trilinos_ENABLE_OpenMP=ON \
      -D TPL_ENABLE_MPI=ON \
      -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -D Kokkos_ENABLE_COMPILER_WARNINGS=ON \
      -D Amesos2_ENABLE_Basker=ON \
      -D Tpetra_INST_SERIAL:BOOL=ON \
      -D Trilinos_ENABLE_TESTS=OFF \
      ../trilinos_source
