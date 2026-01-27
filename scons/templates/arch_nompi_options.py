
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

# No-MPI configuration for Arch Linux.
# Imports the full arch_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   sudo pacman -Sy python python-numpy python-scipy python-sympy python-matplotlib
#   sudo pacman -Sy gcc scons cmake
#   sudo pacman -Sy boost boost-libs suitesparse hdf5 netcdf lapack zlib metis

from arch_options import *

# Disable MPI
mpi = 'none'
mpi4py = False

# Don't build Trilinos (optional, but faster build without MPI)
build_trilinos = 'never'
trilinos = False
