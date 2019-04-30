#!/bin/bash

set -x -e

INCLUDE_PATH="${PREFIX}/include"
LIBRARY_PATH="${PREFIX}/lib"

CMAKE_INSTALL_PREFIX=${PREFIX}
INCLUDE_DIRECTORIES="${PREFIX}/include"
LDFLAGS=-L"${PREFIX}/lib"

if [ "$(uname)" == "Linux" ]; then

    # python -m pip install svn
    # svn co https://svn.geocomp.uq.edu.au/svn/esys13/trunk ${SRC_DIR}/escript

    cd ${SRC_DIR}/escript
    scons -j"${CPU_COUNT}" \
        options_file="${SRC_DIR}/escript/scons/templates/stretch_options.py" \
        prefix="${PREFIX}" \
        cxx_extra="-fPIC" \
        boost_prefix="${PREFIX}" \
        boost_libs='boost_python27' \
        pythoncmd=`which python` \
        pythonlibpath="${PREFIX}/lib" \
        pythonincpath="${PREFIX}/include/python${PY_VER}" \
        pythonlibname="python2.7" \
        trilinos=1 \
        trilinos_prefix="${PREFIX}" \
        umfpack=0 \
        umfpack_prefix=["${PREFIX}/include","${PREFIX}/lib"] \
        lapack=0 \
        lapack_prefix=["${PREFIX}/include/atlas","${PREFIX}/lib"] \
        lapack_libs=['lapack'] \
        netcdf_prefix=["${PREFIX}/include","${PREFIX}/lib"] \
        netcdf_libs=['netcdf_c++4','netcdf'] \
        werror=0

fi