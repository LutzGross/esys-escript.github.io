
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

# No-MPI configuration for OpenSUSE.
# Imports the full opensuse_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   sudo zypper in python3-devel python3-numpy python3-scipy python3-sympy python3-matplotlib
#   sudo zypper in gcc gcc-c++ gcc-fortran scons cmake
#   sudo zypper in libboost_python3-devel libboost_numpy3-devel libboost_random-devel
#   sudo zypper in hdf5-devel netcdf-devel suitesparse-devel lapack-devel zlib-devel metis-devel

from opensuse_options import *

# Disable MPI
mpi = 'none'
mpi4py = False

# Don't build Trilinos (optional, but faster build without MPI)
build_trilinos = 'never'
trilinos = False
