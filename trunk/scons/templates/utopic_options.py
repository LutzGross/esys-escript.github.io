
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

# This is a template configuration file for escript on Ubuntu Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
boost_libs = ['boost_python-py27']
#mpi = 'OPENMPI'
mpi_prefix = '/usr/lib/openmpi'
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
netcdf = True
#umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
dudley_assemble_flags = '-funroll-loops'

