
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# This is a template configuration file for escript on Windows 10.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
netcdf = 4

import os, subprocess, sys

# cc_flags for boost
# /EHsc : for unresolved external symbol "void __cdecl boost::throw_exception(...
# /MD : for fatal error C1189: #error:  "Mixing a dll boost library with a static runtime is a really bad idea..."
# /DBOOST_ALL_NO_LIB : for https://github.com/boostorg/python/issues/193
# /wd4068 : for suppressing "warning C4068: unknown pragma" caused by "#pragma clang ..."
cc_flags = '/EHsc /MD /DBOOST_ALL_NO_LIB /wd4068'
# cxx_extra = ''
omp_flags = '/openmp'
username = os.environ['USERNAME']

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = ''

cppunit_prefix = 'C:/Users/{usr}/vcpkg/packages/cppunit_x64-windows'.format(usr=username)
cppunit_libs = ['cppunit_dll']

conda_prefix = os.environ.get('CONDA_PREFIX')
if os.environ.get('CONDA_BUILD'):
    pythoncmd = os.environ['PYTHON']
    build_dir = os.environ['BUILD_PREFIX']
    conda_prefix = os.environ['PREFIX']
    lib_prefix = os.environ['LIBRARY_PREFIX']
elif conda_prefix:
    pythoncmd = conda_prefix + '\\python.exe'
    lib_prefix = conda_prefix + '\\Library'

if conda_prefix:
    ld_extra = ' '.join(filter(None, ('/LIBPATH:{cp}\\Library\\lib /LIBPATH:"{cp}\\libs"'.format(
        cp=conda_prefix), ld_extra)))
    netcdf_prefix = lib_prefix
    netcdf_libs = ['netcdf-cxx4','netcdf','libhdf5','libcurl','zlib']
    boost_prefix = lib_prefix
    # boost_libs ='boost_python37-vc140-mt-x64-1_67'
#    boost_libs, compression_libs = [], []
    boost_libs = []
    cmd = "import sys;print(''.join(str(i) for i in sys.version_info[:2]))"
    py_ver = subprocess.Popen([pythoncmd, '-c', cmd], stdout=subprocess.PIPE).stdout.read().strip()
    if isinstance(py_ver, bytes):
        py_ver = py_ver.decode(sys.stdout.encoding)
    for l in os.listdir(boost_prefix + '\\lib'):
        if l.startswith('boost_python'+py_ver) or l.startswith('boost_numpy'+py_ver):
            boost_libs.append(os.path.splitext(l)[0])
#        elif l.startswith('boost_iostreams'):
#            compression_libs.append(os.path.splitext(l)[0])
    # list comprehension not working with scons?
    # boost_libs = [os.path.splitext(l)[0] for l in os.listdir(boost_prefix + '\\lib')
    #     if l.startswith('boost_python')][-1]
    mumps_seq = True
    mumps_seq_prefix = lib_prefix + '\\mingw-w64'
    ld_extra = ' '.join(filter(None, ('/LIBPATH:"{mp}\\lib" /LIBPATH:"{mp}\\bin"'.format(mp=mumps_seq_prefix), ld_extra)))
    ld_extra = ' '.join(filter(None, ('libmumps_common.a libdmumps.dll.a libzmumps.dll.a', ld_extra)))
#    mumps_seq_libs = ['libmumps_common.a', 'libdmumps.dll.a', 'libzmumps.dll.a']
    mumps_seq_libs = []
    silo_prefix = lib_prefix
    silo_libs = ['silohdf5']

tools_names = ['msvc']
