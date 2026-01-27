
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

# No-MPI configuration for Fedora Linux.
# Imports the full fedora_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   sudo dnf install python3-devel python3-numpy python3-scipy python3-sympy python3-matplotlib
#   sudo dnf install gcc-c++ gcc-gfortran scons cmake
#   sudo dnf install boost-devel boost-python3-devel boost-python3 boost-numpy3
#   sudo dnf install hdf5-devel netcdf-devel suitesparse-devel lapack-devel zlib-devel metis-devel

from fedora_options import *

# Disable MPI
mpi = 'none'
mpi4py = False

# Don't build Trilinos (optional, but faster build without MPI)
build_trilinos = 'never'
trilinos = False
