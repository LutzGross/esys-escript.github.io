
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

# No-MPI configuration for macOS with MacPorts.
# Imports the full macports_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   sudo port install python311
#   sudo port select --set python python311
#   sudo port select --set python3 python311
#   sudo port install py311-numpy py311-scipy py311-sympy py311-matplotlib
#   sudo port install scons cmake boost hdf5 netcdf suitesparse lapack

from scons.templates.macports_options import *

# Disable MPI
mpi = 'none'
mpi4py = False
