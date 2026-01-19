
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
#cxx_extra = '-Wno-literal-suffix'
openmp = True
# mpi='OPENMPI'

prefix='/app'

import os
import sys

cxx_extra='-w -O3 -march=native'

#pythoncmd='/app/bin/python3'
pythonlibpath = '/app/lib'
pythonincpath = '/usr/include/python3.12'

boost_prefix=['/app/include','/app/lib']
boost_libs = ['boost_python312','boost_iostreams','boost_random']
domains=['finley','ripley','speckley']
hdf5_prefix=['/app/include','/app/lib']
hdf5_libs=['hdf5_cpp','hdf5','hdf5_hl_cpp']
mpi_prefix=['/app/include','/app/lib']
mpi_libs=['mpi_cxx']
lapack=0
lapack_prefix = ['/app/include/', '/app/lib']
ld_extra=''
paso=1
p4est=0
silo=0
silo_prefix = ['/app/include/', '/app/lib']
silo_libs = ['silo','hdf5']
trilinos=True
trilinos_prefix = ['/app/include/', '/app/lib']
trilinos_make_sh='tools/flatpak/flatpak_trilinos.sh'
umfpack=False
umfpack_prefix = ['/app/include','/app/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
visit=0
werror=0


# boost-python library/libraries to link against





