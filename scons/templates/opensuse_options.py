now
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# This is a template configuration file for escript on OpenSUSE.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
boost_libs = ['boost_python-py3']
pythoncmd = '/usr/bin/python3'
pythonlibpath = ['/usr/lib64']
pythonlibname = ['python3.6m']
pythonincpath = ['/usr/include/python3.6m']
mumps_seq_libs=['mumps_common','cmumps_seq','dmumps_seq','zmumps_seq','zmumps_seq','pord']
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib64']
umfpack_libs = ['umfpack', 'amd']