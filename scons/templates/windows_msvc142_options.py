
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
verbose = True

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

# Additional flags to add to the linker
# DEFAULT: '' (empty)
#ld_extra = '/LIBPATH:"C:/Program Files (x86)/Windows Kits/10/Lib/10.0.18362.0/ucrt/x86" \
#            /LIBPATH:"C:/Program Files (x86)/Windows Kits/10/Lib/10.0.18362.0/um/x86"'
ld_extra = '/LIBPATH:"C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/netcdf-c_x64-windows/lib" \
            /LIBPATH:"C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/hdf5_x64-windows/lib" \
            /LIBPATH:"C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/curl_x64-windows/lib" \
            /LIBPATH:"C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/zlib_x64-windows/lib" \
            /LIBPATH:"C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/szip_x64-windows/lib"'

netcdf_prefix = 'C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/netcdf-cxx4_x64-windows'
netcdf_libs = ['netcdf-cxx4','netcdf','libhdf5','libcurl','zlib','szip']

cppunit_prefix = 'C:/Users/Mark/Documents/work/sees/escript/vcpkg/packages/cppunit_x64-windows'
cppunit_libs = ['cppunit_dll']

conda_prefix = os.environ.get('CONDA_PREFIX')
if conda_prefix:
    boost_prefix = conda_prefix + '\\Library'
    # boost_libs ='boost_python37-vc140-mt-x64-1_67'
    boost_libs = []
    for l in os.listdir(boost_prefix + '\\lib'):
        if l.startswith('boost_python{}{}'.format(sys.version_info.major,sys.version_info.minor)):
            boost_libs.append(os.path.splitext(l)[0])
    # list comprehension not working with scons?
    # boost_libs = [os.path.splitext(l)[0] for l in os.listdir(boost_prefix + '\\lib') if l.startswith('boost_python')][-1]
    # netcdf_prefix = conda_prefix + '\\Library'
    # netcdf_libs = ['netcdf']

tools_names = ['msvc']
