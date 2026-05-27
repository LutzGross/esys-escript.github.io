#!/bin/bash

set -x -e
set -o pipefail

CFLAGS="${CFLAGS} -I${PREFIX}/include -fPIC"
CXXFLAGS="${CXXFLAGS} -i sysroot ${CONDA_BUILD_SYSROOT} -fPIC -w -fopenmp"

python --version
find ${PREFIX} -iname libboost_python* 2>/dev/null
find ${PREFIX} -iname Python.h 2>/dev/null
find ${PREFIX} -iname libpython* 2>/dev/null
find ${PREFIX} -iname libumfpack* 2>/dev/null

BOOST_LIBS="boost_python${CONDA_PY}"
PYTHON_LIB_PATH="${PREFIX}/lib"
PYTHON_INC_PATH="${PREFIX}/include/python${PY_VER}"
PYTHON_LIB_NAME="python${PY_VER}"
BUILD_SILO=0

cd ${SRC_DIR}/escript
scons -j"${CPU_COUNT}" \
    options_file="${SRC_DIR}/escript/scons/templates/anaconda_options.py" \
    build_dir=${BUILD_PREFIX}/escript_build \
    boost_prefix=${PREFIX} \
    boost_libs=${BOOST_LIBS} \
    cxx=${CXX} \
    cxx_extra="-w -fPIC -fdiagnostics-color -std=c++17 " \
    cppunit_prefix=${PREFIX} \
    ld_extra="-L${PREFIX}/lib -lgomp" \
    netcdf='no' \
    openmp=1 \
    omp_flags="-fopenmp" \
    prefix=${PREFIX} \
    pythoncmd=${PREFIX}/bin/python \
    pythonlibpath=${PYTHON_LIB_PATH} \
    pythonincpath=${PYTHON_INC_PATH} \
    pythonlibname=${PYTHON_LIB_NAME} \
    silo=${BUILD_SILO} \
    silo_prefix=${PREFIX} \
    umfpack_prefix=${PREFIX} \
    verbose=1 \
    build_full || cat config.log

cp -R ${SRC_DIR}/escript/LICENSE ${SRC_DIR}/LICENSE
cp -R ${PREFIX}/esys ${SP_DIR}/esys
cp -R ${BUILD_PREFIX}/escript_build/scripts/release_sanity.py /tmp/release_sanity.py
