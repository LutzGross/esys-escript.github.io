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
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

boost_libs=['boost_python311','boost_numpy311','boost_iostreams']
cxx_extra=' -fdiagnostics-color=always -std=c++17 -Wno-mismatched-new-delete'
escript_opts_version = 203
ld_extra=' -Wl,-soname=libescript.5.10'
mpi4py=0
mumps=1
mumps_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
mumps_libs=['cmumps','dmumps','smumps','zmumps','mumps_common','pord']
netcdf=4
openmp = True
pythonlibname=['python3.11']
pythonincpath='/usr/include/python3.11/'
umfpack=1
umfpack_prefix=['/usr/include/suitesparse/','/usr/lib']
silo=1
silo_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
silo_libs=['siloh5']
verbose=1
