
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

build_trilinos='make'
escript_opts_version=203
mpi='OPENMPI'
mpi_prefix = ['/usr/include/x86_64-linux-gnu/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi' ]
openmp=1
paso=1
cxx_extra= ['-fdiagnostics-color=always', '-Wno-format-truncation']
pythoncmd = 'python3'

parmetis = True

import subprocess
p=subprocess.run([pythoncmd,'-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1]
pythonlibname = 'python3.%s'%subversion
pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = '/usr/include/python3.%s'%subversion
boost_libs=['boost_python3%s'%subversion,'boost_numpy3%s'%subversion,'boost_random']
boost_libs='boost_python3%s'%subversion
boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
werror=0

umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

silo = True
silo_libs = ['siloh5', 'hdf5']
silo_prefix=[ '/usr/include' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']
netcdf = 4

