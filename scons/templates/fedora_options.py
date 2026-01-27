
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

# This is a template configuration file for escript on Fedora Linux.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   sudo dnf install python3-devel python3-numpy python3-scipy python3-sympy python3-matplotlib
#   sudo dnf install gcc-c++ gcc-gfortran scons cmake
#   sudo dnf install boost-devel boost-python3-devel boost-python3 boost-numpy3
#   sudo dnf install hdf5-devel netcdf-devel suitesparse-devel lapack-devel zlib-devel metis-devel
#   sudo dnf install openmpi-devel python3-mpi4py-openmpi
#
# For no-MPI build, use fedora_nompi_options.py instead.

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

# Boost configuration
boost_libs = [f'boost_python3{subversion}']
boost_prefix = ['/usr/include', '/usr/lib64']
disable_boost_numpy = True  # Fedora boost-numpy may need separate installation

# MPI configuration
# Note: On Fedora, you may need to load the MPI module first:
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

# NetCDF configuration
netcdf = True
netcdf_prefix = ['/usr/include', '/usr/lib64']
netcdf_libs = ['netcdf_c++4', 'netcdf']

# HDF5 configuration
hdf5 = True
hdf5_prefix = ['/usr/include', '/usr/lib64']
hdf5_libs = ['hdf5_cpp', 'hdf5']

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

# Optional features
sympy = True
