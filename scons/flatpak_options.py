
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
#cxx_extra = '-Wno-literal-suffix'
openmp = True
# mpi='OPENMPI'
umfpack=True
silo=True
trilinos=True

prefix='/app'

import os
import sys

cxx_extra='-w -O3'

#pythoncmd='/app/bin/python3'
pythonlibpath = '/app/lib'
pythonincpath = '/use/include/python3.12'

boost_prefix=['/app/include','/app/lib']
hdf5_prefix=['/app/include','/app/lib']
mpi_prefix=['/app/include','/app/lib']
mpi_libs=['mpi_cxx']
umfpack_prefix = ['/app/include','/app/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
lapack_prefix = ['/app/include/', '/app/lib']
silo_prefix = ['/app/include/', '/app/lib']
silo_libs = ['silo','hdf5']
trilinos_prefix = ['/app/include/', '/app/lib']

dudley_assemble_flags = '-funroll-loops'

# boost-python library/libraries to link against
boost_libs = ['boost_python312','boost_iostreams','boost_random','boost_numpy312']

trilinos_make_sh='tools/flatpak/flatpak_trilinos.sh'

paso=1
werror=0
