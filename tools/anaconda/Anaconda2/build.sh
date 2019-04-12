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

CMAKE_INSTALL_PREFIX=${PREFIX}
INCLUDE_DIRECTORIES="${PREFIX}/include"
LDFLAGS=-L"${PREFIX}/lib"

if [ "$(uname)" == "Linux" ]; then

    # None of the gcc libraries on Anaconda contain libm.so
    cd ${SRC_DIR}/openlibm_source
    sed -i 's+prefix = /usr/local+prefix = ${PREFIX}+g' Make.inc
    make -j"${CPU_COUNT}" install | tee openlibm.log 2>&1 
    ln -s ${PREFIX}/lib/libopenlibm.so ${PREFIX}/lib/libm.so
    ln -s ${PREFIX}/lib/libopenlibm.so.1 ${PREFIX}/lib/libm.so.1
    ln -s ${PREFIX}/lib/libopenlibm.so.1.0 ${PREFIX}/lib/libm.so.1.0
    ln -s ${PREFIX}/lib/libopenlibm.a ${PREFIX}/lib/libm.a

    mkdir ${SRC_DIR}/lapack_build
    cd ${SRC_DIR}/lapack_build
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_INSTALL_LIBDIR="${PREFIX}/lib" \
        -DCMAKE_CC_FLAGS='-fPIC' \
        -DCMAKE_CXX_FLAGS='-fPIC' \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_TESTING=OFF \
        -DLAPACKE=ON \
        -DCBLAS=ON \
        ${SRC_DIR}/lapack_source
    make -j"${CPU_COUNT}" install | tee lapack.log 2>&1

    # mkdir ${SRC_DIR}/netcdf_build
    # cd ${SRC_DIR}/netcdf_build
    # cmake \
    #      -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    #      -DCMAKE_SHARED_LINKER_FLAGS=${LDFLAGS} \
    #      -DCMAKE_CXX_FLAGS='-fPIC' \
    #      -DBUILD_SHARED_LIBS=ON \
    #      ${SRC_DIR}/netcdf_source 
    # make -j"${CPU_COUNT}" install | tee netcdf.log 2>&1

    # cd ${SRC_DIR}/netcdf_cxx_source
    # ./configure --prefix=${PREFIX}
    # make -j"${CPU_COUNT}" install | tee netcdf_cxx.log 2>&1

    cd ${SRC_DIR}/boost_source
    ./bootstrap.sh \
        --prefix="${PREFIX}" \
        --with-python="${PYTHON}" \
        --with-python-root="${PREFIX} : "${PREFIX}/include"/python${PY_VER}m "${PREFIX}/include"/python${PY_VER}" \
        --with-icu="${PREFIX}" | tee bootstrap.log 2>&1
    ./b2 -q -d+2 -a \
        variant=release \
        address-model="${ARCH}" \
        architecture=x86 \
        debug-symbols=off \
        threading=multi \
        runtime-link=shared \
        link=shared \
        toolset=gcc \
        cxxflags='-fPIC -w' \
        python="${PY_VER}" \
        include=""${PREFIX}/include"" \
        linkflags="-L"${PREFIX}/lib"" \
        --with-python \
        --with-iostreams \
        --with-random \
        --layout=system \
        -j"${CPU_COUNT}" \
        install | tee b2.log 2>&1

    # trilinos
    cd ${SRC_DIR}/trilinos_source
    git reset --hard c87c9f7224b36c845aa7ef19a636f9c99c9c02d4

    mkdir ${SRC_DIR}/trilinos_build
    cd ${SRC_DIR}/trilinos_build
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_MAKE_PROGRAM=${PREFIX}/bin/make \
        -DCMAKE_VERBOSE_MAKEFILE=FALSE \
        -DTrilinos_ENABLE_CXX11=ON \
        -DTrilinos_ENABLE_Fortran=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DOpenMP_C_FLAGS='-w -fopenmp ' \
        -DOpenMP_CXX_FLAGS='-w -fopenmp ' \
        -DCMAKE_CXX_FLAGS='-fPIC -L${PREFIX}/lib -lblas -llapack -lumfpack' \
        -DCMAKE_CC_FLAGS='-fPIC -L${PREFIX}/lib -lblas -llapack -lumfpack' \
        -DTPL_ENABLE_BLAS=ON \
        -DBLAS_LIBRARY_NAMES='blas' \
        -DTPL_BLAS_INCLUDE_DIRS="${PREFIX}/include" \
        -DTPL_BLAS_LIBRARY_DIRS="${PREFIX}/lib" \
        -DTPL_BLAS_LIBRARY_NAMES=${PREFIX}/lib/libblas.so \
        -DTPL_ENABLE_LAPACK=ON \
        -DTPL_LAPACK_INCLUDE_DIRS="${PREFIX}/include" \
        -DTPL_LAPACK_LIBRARY_DIRS="${PREFIX}/lib" \
        -DTPL_LAPACK_LIBRARY_NAMES=${PREFIX}/lib/liblapack.so \
        -DTPL_ENABLE_Boost=ON \
        -DTPL_Boost_INCLUDE_DIRS="${PREFIX}/include" \
        -DTPL_Boost_LIBRARY_DIRS="${PREFIX}/lib" \
        -DTPL_ENABLE_Cholmod=ON \
        -DTPL_Cholmod_INCLUDE_DIRS="${PREFIX}/include" \
        -DTPL_Cholmod_LIBRARY_DIRS="${PREFIX}/lib" \
        -DTPL_Cholmod_LIBRARIES='libcholmod.so;libamd.so;libcolamd.so' \
        -DTPL_ENABLE_METIS=ON \
        -DTPL_METIS_INCLUDE_DIRS="${PREFIX}/include" \
        -DTPL_METIS_LIBRARY_DIRS="${PREFIX}/lib" \
        -DTPL_ENABLE_UMFPACK=ON \
        -DTPL_UMFPACK_INCLUDE_DIRS="${PREFIX}/include" \
        -DTPL_UMFPACK_LIBRARY_DIRS="${PREFIX}/lib" \
        -DTrilinos_ENABLE_Amesos2=ON \
        -DTrilinos_ENABLE_Belos=ON \
        -DTrilinos_ENABLE_Ifpack2=ON \
        -DTrilinos_ENABLE_Kokkos=ON \
        -DTrilinos_ENABLE_MueLu=ON \
        -DTrilinos_ENABLE_Tpetra=ON \
        -DTrilinos_ENABLE_Teuchos=ON \
        -DTrilinos_ENABLE_COMPLEX=ON \
        -DTrilinos_ENABLE_OpenMP=ON \
        -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
        -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
        ${SRC_DIR}/trilinos_source
    make -j"${CPU_COUNT}" install

    cd ${SRC_DIR}/escript
    scons -j"${CPU_COUNT}" \
        options_file="${SRC_DIR}/escript/scons/templates/stretch_options.py" \
        prefix="${PREFIX}" \
        cxx_extra="-fPIC" \
        ld_extra="-L${PREFIX}/lib -lblas -llapack -lumfpack" \
        boost_prefix="${PREFIX}" \
        boost_libs='boost_python27' \
        pythoncmd=`which python` \
        pythonlibpath="${PREFIX}/lib" \
        pythonincpath="${PREFIX}/include/python${PY_VER}" \
        pythonlibname="python2.7" \
        trilinos=1 \
        trilinos_prefix="${PREFIX}" \
        umfpack=0 \
        umfpack_prefix=[""${PREFIX}/include"",""${PREFIX}/lib""] \
        lapack=0 \
        lapack_prefix=[""${PREFIX}/include"/atlas",""${PREFIX}/lib""] \
        lapack_libs=['lapack'] \
        werror=0

fi