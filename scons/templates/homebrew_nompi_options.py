
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

# No-MPI configuration for macOS with Homebrew.
# Imports the full homebrew_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   brew install python3 scons cmake llvm
#   brew install boost boost-python3
#   brew install hdf5 suite-sparse netcdf netcdf-cxx4
#   pip3 install numpy scipy sympy matplotlib

from homebrew_options import *

# Disable MPI
mpi = 'none'
mpi4py = False

# Don't build Trilinos (optional, but faster build without MPI)
build_trilinos = 'never'
trilinos = False
