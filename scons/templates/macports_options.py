
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

# This is a template configuration file for escript on macOS with MacPorts.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   sudo port install python311
#   sudo port select --set python python311
#   sudo port select --set python3 python311
#   sudo port install py311-numpy py311-scipy py311-sympy py311-matplotlib
#   sudo port install scons cmake boost hdf5 netcdf suitesparse lapack
#   sudo port install openmpi py311-mpi4py
#
# For no-MPI build, use macports_nompi_options.py instead.

import subprocess

escript_opts_version = 203

# MacPorts prefix
MACPORTS_PREFIX = '/opt/local'

# Python configuration - auto-detect version
pythoncmd = '/opt/local/bin/python3'
p = subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
subversion = p.stdout.split(' ')[1].split('.')[1]
print(f"Python subversion = {subversion}")

pythonlibpath = [f'{MACPORTS_PREFIX}/Library/Frameworks/Python.framework/Versions/3.{subversion}/lib']
pythonincpath = [f'{MACPORTS_PREFIX}/Library/Frameworks/Python.framework/Versions/3.{subversion}/include/python3.{subversion}']

# Compiler settings
openmp = True
tools_names = ['clang']

# MPI configuration
mpi = 'OPENMPI'
mpi_prefix = [f'{MACPORTS_PREFIX}/include/openmpi-mp', f'{MACPORTS_PREFIX}/lib/openmpi-mp']
mpi_libs = ['mpi_cxx', 'mpi']
mpi4py = True

# Solver configuration
paso = True
build_trilinos = 'check'

# Boost configuration
boost_prefix = MACPORTS_PREFIX
compression_libs = ['boost_iostreams-mt']

# NetCDF configuration
netcdf = True
netcdf_prefix = MACPORTS_PREFIX

# SILO configuration
silo = True
silo_prefix = MACPORTS_PREFIX
silo_libs = ['siloh5', 'hdf5']

# HDF5 configuration
hdf5 = True
hdf5_prefix = MACPORTS_PREFIX
hdf5_libs = ['hdf5_cpp', 'hdf5']

# UMFPACK direct solver
umfpack = True
umfpack_prefix = [f'{MACPORTS_PREFIX}/include/suitesparse', f'{MACPORTS_PREFIX}/lib']

# LAPACK configuration
lapack = 'auto'

# zlib configuration
zlib = True
zlib_prefix = MACPORTS_PREFIX
zlib_libs = ['z']

# Optional features
sympy = True

# Clean up
del MACPORTS_PREFIX
