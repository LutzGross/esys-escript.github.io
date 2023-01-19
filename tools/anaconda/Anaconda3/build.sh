#!/bin/bash

set -x -e
set -o pipefail

CXXFLAGS="${CXXFLAGS} -fPIC -w -fopenmp"

cd ${SRC_DIR}/escript
if [ ${PY3K} -eq 1 ]
then
    scons -j"${CPU_COUNT}" \
        options_file="${SRC_DIR}/escript/scons/templates/anaconda_python3_options.py" \
        build_dir=${BUILD_PREFIX}/escript_build \
        boost_prefix=${PREFIX} \
        boost_libs="boost_python${CONDA_PY}" \
        cxx=${CXX} \
        cxx_extra="-w -fPIC" \
        cppunit_prefix=${PREFIX} \
        ld_extra="-L${PREFIX}/lib -lgomp" \
        omp_flags="-fopenmp" \
        prefix=${PREFIX} \
        pythoncmd=${PREFIX}/bin/python \
        pythonlibpath=${PREFIX}/lib \
        pythonincpath=${PREFIX}/include/python${PY_VER} \
        pythonlibname=python${PY_VER} \
        umfpack_prefix=${PREFIX} \
        build_full || cat config.log
else
    scons -j"${CPU_COUNT}" \
        options_file="${SRC_DIR}/escript/scons/templates/buster_py2_options.py" \
        prefix=${PREFIX} \
        build_dir=${BUILD_PREFIX}/escript_build \
        boost_prefix=${PREFIX} \
        boost_libs="boost_python27" \
        cxx=${CXX} \
        cxx_extra="-w -fPIC" \
        ld_extra="-L${PREFIX}/lib -lgomp" \
        cppunit_prefix=${PREFIX} \
        openmp=1 \
        omp_flags="-fopenmp" \
        pythoncmd=${PREFIX}/bin/python \
        pythonlibpath="${PREFIX}/lib" \
        pythonincpath="${PREFIX}/include/python2.7" \
        pythonlibname="python2.7" \
        paso=1 \
        trilinos=0 \
        umfpack=0 \
        umfpack_prefix="${PREFIX}" \
        netcdf=no \
        werror=0 \
        verbose=0 \
        compressed_files=0 \
        build_full || cat config.log
fi

cp -R ${SRC_DIR}/escript/LICENSE ${SRC_DIR}/LICENSE
cp -R ${PREFIX}/esys ${SP_DIR}/esys
cp -R ${BUILD_PREFIX}/escript_build/scripts/release_sanity.py /tmp/release_sanity.py
