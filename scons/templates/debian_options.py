
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

# This is a template configuration file for escript on Debian/Ubuntu Linux.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   sudo apt-get install python3-dev python3-numpy python3-scipy python3-sympy
#   sudo apt-get install scons cmake g++ gfortran
#   sudo apt-get install libboost-python-dev libboost-numpy-dev libboost-random-dev libboost-iostreams-dev
#   sudo apt-get install libnetcdf-dev libsilo-dev libhdf5-serial-dev
#   sudo apt-get install libsuitesparse-dev liblapacke-dev libmumps-seq-dev
#   sudo apt-get install libmetis-dev libparmetis-dev zlib1g-dev
#   sudo apt-get install libopenmpi-dev openmpi-bin python3-mpi4py
#
# For no-MPI build, use debian_nompi_options.py instead.

import subprocess

escript_opts_version = 203

# Python configuration - auto-detect version
pythoncmd = 'python3'
p = subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
subversion = p.stdout.split(' ')[1].split('.')[1]
print(f"Python subversion = {subversion}")

pythonlibname = f'python3.{subversion}'
pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = f'/usr/include/python3.{subversion}'

# Compiler settings
openmp = True
werror = False
cxx_extra = ['-O3', '-fdiagnostics-color=always', '-fstack-protector-strong', '-Wformat', '-Werror=format-security']

# Boost configuration
boost_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu/']
boost_libs = [f'boost_python3{subversion}', f'boost_numpy3{subversion}', 'boost_random']

# MPI configuration
# Note: OpenMPI 5.x removed mpi_cxx and open-rte libraries; just 'mpi' is needed
mpi = 'OPENMPI'
mpi_prefix = ['/usr/include/x86_64-linux-gnu/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi']
mpi_libs = ['mpi']
mpi4py = True

# Solver configuration
paso = True
build_trilinos = 'check'

# UMFPACK direct solver
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

# MUMPS sequential solver
mumps_seq = True
mumps_seq_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
mumps_seq_libs = ['dmumps_seq', 'zmumps_seq', 'mumps_common_seq', 'mpiseq_seq', 'pord_seq']

# METIS graph partitioning (for Trilinos)
metis = True
metis_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
metis_libs = ['metis']

# ParMETIS parallel graph partitioning (requires MPI)
parmetis = True
parmetis_prefix = ['/usr/include/parmetis', '/usr/lib/x86_64-linux-gnu']
parmetis_libs = ['parmetis']

# LAPACK configuration
lapack = 'auto'
lapack_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
lapack_libs = ['lapacke']

# File I/O libraries
netcdf = True
netcdf_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
netcdf_libs = ['netcdf_c++4', 'netcdf']

silo = True
silo_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']
silo_libs = ['siloh5', 'hdf5']

hdf5 = True
hdf5_prefix = ['/usr/include/hdf5/serial', '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']
hdf5_libs = ['hdf5_serial_cpp']

# Compression support
zlib = True
zlib_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
zlib_libs = ['z']

# Optional features
sympy = True
