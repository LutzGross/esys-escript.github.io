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

import sys

py_vers_string = "%s.%s"%(sys.version_info[0],sys.version_info[1])

# locations of include files for python
python_path = "/usr/include/python" + py_vers_string
python_lib_path = '/usr/lib'
python_lib = "python" + py_vers_string

# locations of libraries for boost
boost_path = '/usr/include'
boost_lib_path = '/usr/lib'
boost_lib = 'boost_python'

useNetCDF = 'yes'
netCDF_path = '/usr/include'
netCDF_lib_path = '/usr/lib'
netCDF_libs = ['netcdf', 'netcdf_c++']

# locations of doc building executables
doxygen_path = '/usr/bin'
epydoc_path = '/usr/bin'

# c flags to use
cc_flags  = '-fpic  -Wno-unknown-pragmas'
cc_optim  = '-O3'
cc_debug  = '-g'

# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
#sys_libs = []
