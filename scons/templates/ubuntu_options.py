
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
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

build_trilinos='make'
escript_opts_version=203
mpi='OPENMPI'
mpi_prefix = ['/usr/include/x86_64-linux-gnu/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi' ]
openmp = 1
paso = 1
sympy = True
#cxx_extra = ['-freference fixed. diagnostics-color=always', '-Wno-format-truncation']
cxx_extra= ['-O3', '-fdiagnostics-color=always', '-fstack-protector-strong',  '-Wformat', '-Werror=format-security' ]
pythoncmd = 'python3'

parmetis = True

import subprocess
p=subprocess.run([pythoncmd,'-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1]
pythonlibname = 'python3.%s'%subversion
pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = '/usr/include/python3.%s'%subversion
#boost_libs=['boost_python3%s'%subversion,'boost_numpy3%s'%subversion,'boost_random']
boost_libs=['boost_python3%s'%subversion, 'boost_random']
boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
werror=0

umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

silo = True
silo_libs = ['siloh5', 'hdf5']
silo_prefix=[ '/usr/include' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']


hdf5 = True
hdf5_libs = ['hdf5_serial_cpp']
hdf5_prefix=[ '/usr/include/hdf5/serial' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']

# LAPACK configuration - uses LAPACKE (modern C interface)
# Install with: sudo apt-get install liblapacke-dev
lapack = 'auto'  # Auto-detect LAPACKE
lapack_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
lapack_libs = ['lapacke']

# MUMPS configuration - sequential version (works with MPI builds)
# Install with: sudo apt-get install libmumps-seq-dev
mumps_seq = True
mumps_seq_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
mumps_seq_libs = ['dmumps_seq', 'zmumps_seq', 'mumps_common_seq', 'mpiseq_seq', 'pord_seq']
# NetCDF configuration - version 4
# Install with: sudo apt-get install libnetcdf-dev libnetcdf-c++4-dev
netcdf = True
netcdf_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
netcdf_libs = ['netcdf_c++4', 'netcdf']
