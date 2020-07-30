
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

# This is a template configuration file for escript on Windows 10.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
netcdf = 4

import os
import sys

# cc_flags for boost
# /EHsc : for unresolved external symbol "void __cdecl boost::throw_exception(...
# /MD : for fatal error C1189: #error:  "Mixing a dll boost library with a static runtime is a really bad idea..."
# /DBOOST_ALL_NO_LIB : for https://github.com/boostorg/python/issues/193
# /wd4068 : for suppressing "warning C4068: unknown pragma" caused by "#pragma clang ..."
cc_flags = '/EHsc /MD /DBOOST_ALL_NO_LIB /wd4068'
# cxx_extra = ''
omp_flags = '/openmp'
username = os.environ['USERNAME']
vcpkg_prefix = 'C:/Users/{usr}/vcpkg/packages'.format(usr=username)

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = '/LIBPATH:"{vcp}/hdf5_x64-windows/lib" /LIBPATH:"{vcp}/curl_x64-windows/lib" \
            /LIBPATH:"{vcp}/zlib_x64-windows/lib" /LIBPATH:"{vcp}/szip_x64-windows/lib"'.format(vcp=vcpkg_prefix)

netcdf_prefix = '{vcp}/netcdf-cxx4_x64-windows'.format(vcp=vcpkg_prefix)
netcdf_libs = ['netcdf-cxx4','netcdf','libhdf5','libcurl','zlib','szip']
ld_extra = ' /LIBPATH:"{ncp}/lib" {ld}'.format(ncp=netcdf_prefix, ld=ld_extra)

cppunit_prefix = '{vcp}/cppunit_x64-windows'.format(vcp=vcpkg_prefix)
cppunit_libs = ['cppunit_dll']

conda_prefix = os.environ['CONDA_PREFIX']
if conda_prefix:
    ld_extra = ' /LIBPATH:"{cp}\\libs" {ld}'.format(cp=conda_prefix, ld=ld_extra)
    boost_prefix = conda_prefix + '\\Library'
    # boost_libs ='boost_python37-vc140-mt-x64-1_67'
    boost_libs = []
    for l in os.listdir(boost_prefix + '\\lib'):
        if l.startswith('boost_python{}{}'.format(sys.version_info.major,sys.version_info.minor)):
            boost_libs.append(os.path.splitext(l)[0])
    # list comprehension not working with scons?
    # boost_libs = [os.path.splitext(l)[0] for l in os.listdir(boost_prefix + '\\lib')
    #     if l.startswith('boost_python')][-1]
    mumps = True
    mumps_prefix = conda_prefix + '\\Library\\mingw-w64'
    ld_extra = ' {ld} /LIBPATH:"{mp}\\lib" /LIBPATH:"{mp}\\bin"'.format(mp=mumps_prefix, ld=ld_extra)
    ld_extra = ld_extra+' libmumps_common.a libdmumps.dll.a libzmumps.dll.a'
    mumps_libs = []

tools_names = ['msvc']
