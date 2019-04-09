#!/bin/bash

# Hints:
# http://boost.2283326.n4.nabble.com/how-to-build-boost-with-bzip2-in-non-standard-location-td2661155.html
# http://www.gentoo.org/proj/en/base/amd64/howtos/?part=1&chap=3
# http://www.boost.org/doc/libs/1_55_0/doc/html/bbv2/reference.html

# Hints for OSX:
# http://stackoverflow.com/questions/20108407/how-do-i-compile-boost-for-os-x-64b-platforms-with-stdlibc

set -x -e

INCLUDE_PATH="${PREFIX}/include"
LIBRARY_PATH="${PREFIX}/lib"

if [ "$(uname)" == "Linux" ]; then

    cd ${SRC_DIR}/Boost_Source

    ./bootstrap.sh \
        --prefix="${PREFIX}" \
        --with-python="${PYTHON}" \
        --with-python-root="${PREFIX} : ${PREFIX}/include/python${PY_VER}m ${PREFIX}/include/python${PY_VER}" \
        --with-icu="${PREFIX}" \
        | tee bootstrap.log 2>&1

    ./b2 -q -d+2 -a \
        variant=release \
        address-model="${ARCH}" \
        architecture=x86 \
        debug-symbols=off \
        threading=multi \
        runtime-link=shared \
        link=static,shared \
        toolset=gcc \
        python="${PY_VER}" \
        include="${INCLUDE_PATH}" \
        linkflags="-L${LIBRARY_PATH}" \
        --with-python \
        --with-iostreams \
        --with-random \
        --layout=system \
        -j"${CPU_COUNT}" \
        install | tee b2.log 2>&1

    cd ${SRC_DIR}/Trilinos_Source
    git reset --hard c87c9f7224b36c845aa7ef19a636f9c99c9c02d4

    mkdir ${SRC_DIR}/trilinos_build
    cd ${SRC_DIR}/trilinos_build
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_MAKE_PROGRAM=${PREFIX}/bin/make \
        -DTrilinos_ENABLE_CXX11=ON \
        -DTrilinos_ENABLE_Fortran=OFF \
        -DBUILD_SHARED_LIBS=ON \
        -DOpenMP_C_FLAGS='-w -fopenmp' \
        -DOpenMP_CXX_FLAGS='-w -fopenmp' \
        -DTPL_ENABLE_BLAS=ON \
        -DBLAS_LIBRARY_NAMES='blas' \
        -DTPL_ENABLE_LAPACK=ON \
        -DTPL_ENABLE_Boost=ON \
        -DTPL_ENABLE_Cholmod=OFF \
        -DTPL_ENABLE_Matio=OFF \
        -DTPL_ENABLE_METIS=OFF \
        -DTPL_ENABLE_ParMETIS=OFF \
        -DMETIS_LIBRARY_NAMES='metis' \
        -DTPL_ENABLE_ScaLAPACK=OFF \
        -DTPL_ENABLE_UMFPACK=OFF \
        -DTrilinos_ENABLE_Amesos2=ON \
        -DTrilinos_ENABLE_Belos=ON \
        -DTrilinos_ENABLE_Ifpack2=ON \
        -DTrilinos_ENABLE_Kokkos=ON \
        -DTrilinos_ENABLE_MueLu=ON \
        -DTrilinos_ENABLE_Tpetra=ON \
        -DTrilinos_ENABLE_Teuchos=ON \
        -DTrilinos_ENABLE_COMPLEX=ON \
        -DTrilinos_ENABLE_OpenMP=ON \
        -DTrilinos_ENABLE_MPI=OFF \
        -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON \
        -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
        ${SRC_DIR}/Trilinos_Source
    make -j"${CPU_COUNT}" install

    cd ${SRC_DIR}/escript
    scons -j"${CPU_COUNT}" \
        options_file="${SRC_DIR}/escript/scons/templates/sid_options.py" \
        prefix="${PREFIX}" \
        trilinos=True \
        trilinos_prefix="${PREFIX}" \
        pythonlibpath="${PREFIX}/lib" \
        pythonincpath="${PREFIX}/include/python${PY_VER}" \
        pythonlibname="python2.7" \
        werror=0

fi