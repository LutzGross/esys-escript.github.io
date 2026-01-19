ne
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file in the repository root for contributors and development history
# https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS
#
##############################################################################

# This is a template configuration file for escript on CentOS.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
netcdf = True
compressed_files = False
boost_libs = ['boost_python3']
boost_prefix=['/usr/include/','/usr/lib64']
netcdf='no'
pythoncmd='/usr/bin/python3'
pythonlibpath = ['/usr/lib64']
pythonlibname = ['python3.6m']
pythonincpath = ['/usr/include/python3.6m']
cxx_extra = '-std=c++11'
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse','/usr/lib']

