
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

boost_prefix='/usr/local/boost.1.74.0'
boost_libs='boost_python38'
# cxx="/usr/bin/clang++"
# cxx_extra="-O3 -funroll-loops -fopenmp -fdiagnostics-color"
cxx_extra="-O0 -g -pg -fdiagnostics-color -Wno-implicit-int-float-conversion"
debug=True
# debug=False
# # lapack=True
# lapack_prefix=['/usr/include/x86_64-linux-gnu/','/usr/lib']
ld_extra='-L/usr/lib/openmpi/'
mpi = 'OPENMPI'
mpi_prefix=['/usr/include/','/usr/lib/openmpi/']
mpi_libs=['mpi','mpi_cxx']
# netcdf="no"
# # netcdf=4
openmp=True
# openmp=False
# paso=True
paso=False
pythoncmd="/usr/bin/python3"
pythonlibname="python3.8"
pythonlibpath="/usr/lib"
pythonincpath="/usr/include/python3.8"
umfpack=True
umfpack_prefix=['/usr/include/','/usr/lib']
# # silo=True
# silo_libs=['siloh5']
# silo_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
trilinos=True
# trilinos_prefix='/usr/local/trilinos_noomp'
# trilinos_prefix='/usr/local/trilinos_nompi'
trilinos_prefix='/usr/local/trilinos_mpi'
verbose=True
# visit=True
# visit_libs='simV2runtime_par'
# visit_prefix='/usr/local/visit/3.1.2/linux-x86_64/libsim/V2/'
werror=False