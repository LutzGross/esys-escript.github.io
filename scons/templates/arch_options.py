
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

# This is a template configuration file for escript on Arch Linux.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   sudo pacman -Sy python python-numpy python-scipy python-sympy python-matplotlib
#   sudo pacman -Sy gcc scons cmake
#   sudo pacman -Sy boost boost-libs suitesparse hdf5 netcdf lapack zlib metis
#   sudo pacman -Sy openmpi python-mpi4py
#
# For no-MPI build, use arch_nompi_options.py instead.

import subprocess

escript_opts_version = 203

# Python configuration - auto-detect version
pythoncmd = '/usr/bin/python3'
p = subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
subversion = p.stdout.split(' ')[1].split('.')[1]
print(f"Python subversion = {subversion}")

pythonlibpath = ['/usr/lib']
pythonlibname = [f'python3.{subversion}']
pythonincpath = [f'/usr/include/python3.{subversion}']

# Compiler settings
openmp = True
werror = 0

# Boost configuration
boost_libs = [f'boost_python3{subversion}']
boost_prefix = ['/usr/include', '/usr/lib']
disable_boost_numpy = False

# MPI configuration
mpi = 'OPENMPI'
mpi_prefix = ['/usr/include/openmpi', '/usr/lib/openmpi']
mpi_libs = ['mpi_cxx', 'mpi']
mpi4py = True

# Solver configuration
paso = True
build_trilinos = 'check'

# UMFPACK direct solver
umfpack = True
umfpack_prefix = ['/usr/include/', '/usr/lib']

# HDF5 configuration
hdf5 = True
hdf5_prefix = ['/usr/include', '/usr/lib']
hdf5_libs = ['hdf5_cpp', 'hdf5']

# NetCDF configuration
netcdf = True
netcdf_prefix = ['/usr/include', '/usr/lib']
netcdf_libs = ['netcdf_c++4', 'netcdf']

# LAPACK configuration
lapack = 'auto'

# METIS configuration
metis = True
metis_prefix = ['/usr/include', '/usr/lib']
metis_libs = ['metis']

# zlib configuration
zlib = True
zlib_prefix = ['/usr/include', '/usr/lib']
zlib_libs = ['z']

# Optional features
sympy = True
