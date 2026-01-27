
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

# This is a template configuration file for escript on OpenSUSE.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   sudo zypper in python3-devel python3-numpy python3-scipy python3-matplotlib
#   sudo zypper in gcc gcc-c++ gcc-fortran scons cmake
#   sudo zypper in libboost_python3-devel libboost_numpy3-devel libboost_random-devel
#   sudo zypper in hdf5-devel netcdf-devel suitesparse-devel lapack-devel zlib-devel metis-devel

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

# Boost configuration - OpenSUSE uses different library naming
boost_libs = ['boost_python-py3']
boost_prefix = ['/usr/include', '/usr/lib64']

# Compiler settings
openmp = True
werror = 0

# UMFPACK direct solver
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib64']
umfpack_libs = ['umfpack', 'amd']

# MUMPS sequential solver
mumps_seq = True
mumps_seq_libs = ['mumps_common', 'cmumps_seq', 'dmumps_seq', 'zmumps_seq', 'pord']

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

# MPI disabled by default - enable for parallel support
# sudo zypper in openmpi-devel python3-mpi4py
mpi = 'none'

# Trilinos - build from bundled source
build_trilinos = 'make'
