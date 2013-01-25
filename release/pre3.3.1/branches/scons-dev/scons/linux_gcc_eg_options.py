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

# locations of include files for python
python_path = "/usr/include/python%s.%s"%(sys.version_info[0],sys.version_info[1])
python_lib_path = '/usr/lib'
python_lib = "python%s.%s"%(sys.version_info[0],sys.version_info[1])

# locations of libraries for boost
boost_path = '/usr/include'
boost_lib_path = '/usr/lib'
boost_lib = 'boost_python'

useNetCDF = 'yes'
netCDF_libs = ['netcdf', 'netcdf_c++']

# locations of doc building executables
doxygen_path = '/usr/bin'
epydoc_path = '/usr/bin'

# c flags to use
cc_flags  = '-O3 -fpic  -Wno-unknown-pragmas'
cc_flags_debug  = '-g -O0 -fpic -Wno-unknown-pragmas'

# c++ flags to use
cxx_flags  = '-ansi'
cxx_flags_debug  = '-ansi -DDOASSERT -DDOPROF'

# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
#sys_libs = []
