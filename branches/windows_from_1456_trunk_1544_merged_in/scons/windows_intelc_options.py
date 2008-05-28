#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

import os

source_root = os.path.realpath('.')

pyinstall = os.path.join(source_root,'esys')
incinstall = os.path.join(source_root,'include')
libinstall = os.path.join(source_root,'lib')
exinstall = os.path.join(source_root,'examples')

# locations of include files for python
python_root = 'C:/python25'
python_cmd = os.path.join(python_root,'python')
python_path =  os.path.join(python_root,'include')
python_lib_path = os.path.join(python_root,'libs')
python_lib = 'python25'

# locations of libraries for boost
dotdot = os.path.realpath('..')
boost_path = os.path.join(dotdot,'boost-1.33')
boost_lib_path = os.path.join(boost_path,'windows_binary','lib')
boost_lib = 'boost_python-vc71-mt-s-1_33_1.lib'

# locations of netcdf
useNetCDF = "no"
netCDF_path = os.path.realpath(".")
netCDF_lib_path = os.path.realpath(".")
netCDF_libs = [ ]

cc_defines = ['_USE_MATH_DEFINES', 'DLL_NETCDF', 'BOOST_NO_INTRINSIC_WCHAR_T' ]
# c flags to use
# 1563 - taking adress of a temporary
# 811 - exception specification for implicitly declared virtual function (destructor usually) incompatible with that of override
# 161 - openmp pargmas are unknown when not compiling with openmp
cc_common_flags = '/FD /EHsc /GR /Qvc7.1 '
cc_flags  = cc_common_flags + '/O2 /Op /MT /W3'

cc_flags_debug  = cc_common_flags + '/Od /RTC1 /MTd /ZI'

# c++ flags to use
cxx_flags = ''
cxx_flags_debug = ''
# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
sys_libs = []
