
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


import sys, os

source_root = os.path.realpath('.')

pyinstall = os.path.join(source_root,'esys')
incinstall = os.path.join(source_root,'include')
libinstall = os.path.join(source_root,'lib')
exinstall = os.path.join(source_root,'examples')

# locations of files for python
py_vers = '%s%s'%(sys.version_info[0],sys.version_info[1])
python_root = 'C:/python' + py_vers
python_cmd = os.path.join(python_root,'python')
python_path =  os.path.join(python_root,'include')
python_lib_path = os.path.join(python_root,'libs')
python_libs = ['python' + py_vers]

# locations of libraries for boost
dotdot = os.path.realpath('..')
boost_path = os.path.join(dotdot,'boost-1.33')
boost_lib_path = os.path.join(boost_path,'windows_binary','lib')
boost_libs = ['boost_python-vc71-mt-1_33_1.lib']

# locations of netcdf
useNetCDF = "yes"
netCDF_root = os.path.join(dotdot,"netcdf")
netCDF_path = os.path.join(netCDF_root,"src","include")
netCDF_lib_path = os.path.join(netCDF_root,'lib')
netCDF_libs = ["netcdf", "netcdf_cpp" ]

cc_defines = ['_USE_MATH_DEFINES']
# c flags to use
# /wd4068 - ignore unknown pragmas.
# /EHsc - /EHsc Catch C++ exceptions only, assume C code does not throw.
# /GR - enable RTTI.
# /MD - link to DLL multithreaded RTL
# /ZI - dont put the library name in the .obj file.
# /Yd - Put debug info in .obj file.
# /Y- - Ignore precompiled header options, we don't want them.
# /RTCcsu - Stack frame, use of uninitialised variabe and truncation checks.
cc_flags  = '/EHsc /GR /wd4068  /MD'
cc_optim  = '/O2 /Op /W3'
cc_debug  = '/Od /RTCcsu /ZI /Yd /Y-'

# MPI version using MPICH2
usempi = 'yes'
mpi_path = 'C:/Program Files/MPICH2/include'
mpi_lib_path = 'C:/Program Files/MPICH2/lib'
mpi_libs = ['mpi']
mpi_flavour = 'MPICH2'

# linker flags to use
#link_flags = ''
#link_flags_debug = '/debug /incremental:no /opt:ref /opt:icf'

# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
sys_libs = ["C:/Program Files/Microsoft Visual Studio .NET 2003/Vc7/PlatformSDK/Lib/Ws2_32"]

tools_names = ['msvc']
