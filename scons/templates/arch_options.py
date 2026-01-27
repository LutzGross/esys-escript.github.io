
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
#   sudo pacman -Sy python python-numpy python-scipy python-matplotlib
#   sudo pacman -Sy gcc scons cmake
#   sudo pacman -Sy boost boost-libs suitesparse hdf5 netcdf lapack zlib metis

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

# Boost configuration
boost_libs = [f'boost_python3{subversion}']
boost_prefix = ['/usr/include', '/usr/lib']
disable_boost_numpy = False

# Compiler settings
openmp = True
werror = 0

# UMFPACK direct solver
umfpack = True
umfpack_prefix = ['/usr/include/', '/usr/lib']

# HDF5 configuration
hdf5 = True
hdf5_prefix = ['/usr/include', '/usr/lib']
hdf5_libs = ['hdf5_cpp', 'hdf5']

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

# MPI disabled by default - enable for parallel support
# sudo pacman -Sy openmpi python-mpi4py
mpi = 'none'

# Trilinos - build from bundled source
build_trilinos = 'make'
