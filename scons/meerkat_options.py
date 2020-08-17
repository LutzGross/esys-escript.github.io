
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203

boost_prefix='/usr/local/boost.1.73.0'
boost_libs='boost_python37'
debug=False
lapack=True
lapack_prefix=['/usr/include/x86_64-linux-gnu/','/usr/lib']
ld_extra='-Lhdf5/serial/libhdf5.so'
netcdf=4
openmp=True
pythoncmd='/usr/bin/python3'
pythonlibname='python3.7m'
pythonlibpath='/usr/lib'
pythonincpath='/usr/include/python3.7'
umfpack=True
umfpack_prefix=['/usr/include/suitesparse','/usr/lib']
silo=True
silo_libs=['siloh5']
silo_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
trilinos=True
trilinos_prefix='/usr/local/trilinos_nompi'
# verbose=True
# visit=True
visit_libs='simV2runtime_par'
visit_prefix='/usr/local/visit/3.1.2/linux-x86_64/libsim/V2/'