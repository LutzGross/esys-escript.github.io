
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

# This is a template configuration file for escript on RHEL-compatible distributions:
# RHEL 8/9, Rocky Linux, AlmaLinux, CentOS Stream.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   sudo dnf install epel-release
#   sudo dnf config-manager --set-enabled powertools  # RHEL 8 / Rocky 8
#   # or: sudo dnf config-manager --set-enabled crb   # RHEL 9 / Rocky 9
#   sudo dnf install python3-devel python3-numpy python3-scipy python3-sympy python3-matplotlib
#   sudo dnf install gcc gcc-c++ gcc-gfortran scons cmake
#   sudo dnf install boost-devel boost-python3 boost-python3-devel
#   sudo dnf install hdf5-devel netcdf-devel suitesparse suitesparse-devel lapack-devel zlib-devel metis-devel
#   sudo dnf install openmpi-devel python3-mpi4py-openmpi
#
# For no-MPI build, use rhel_nompi_options.py instead.

import subprocess

escript_opts_version = 203

# Python configuration - auto-detect version
pythoncmd = '/usr/bin/python3'
p = subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
subversion = p.stdout.split(' ')[1].split('.')[1]
print(f"Python subversion = {subversion}")

pythonlibpath = ['/usr/lib64']
pythonlibname = [f'python3.{subversion}']
pythonincpath = [f'/usr/include/python3.{subversion}']

# Compiler settings
openmp = True
werror = 0
cxx_extra = '-std=c++11'

# Boost configuration
boost_libs = ['boost_python3']
boost_prefix = ['/usr/include/', '/usr/lib64']

# MPI configuration
# Note: You may need to load the MPI module first:
#   module load mpi/openmpi-x86_64
mpi = 'OPENMPI'
mpi_prefix = ['/usr/include/openmpi-x86_64', '/usr/lib64/openmpi/lib']
mpi_libs = ['mpi_cxx', 'mpi']
mpi4py = True

# Solver configuration
paso = True
build_trilinos = 'make'

# UMFPACK direct solver
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib64']

# HDF5 configuration
hdf5 = True
hdf5_prefix = ['/usr/include', '/usr/lib64']
hdf5_libs = ['hdf5_cpp', 'hdf5']

# NetCDF - may need manual configuration on CentOS
netcdf = False

# LAPACK configuration
lapack = 'auto'

# METIS configuration
metis = True
metis_prefix = ['/usr/include', '/usr/lib64']
metis_libs = ['metis']

# zlib configuration
zlib = True
zlib_prefix = ['/usr/include', '/usr/lib64']
zlib_libs = ['z']

# Compression support may need boost_iostreams
compressed_files = False

# Optional features
sympy = True
