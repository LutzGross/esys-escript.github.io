
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

# No-MPI configuration for RHEL-compatible distributions:
# RHEL 8/9, Rocky Linux, AlmaLinux, CentOS Stream.
# Imports the full rhel_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   sudo dnf install epel-release
#   sudo dnf config-manager --set-enabled powertools
#   sudo dnf install python3-devel python3-numpy python3-scipy python3-sympy python3-matplotlib
#   sudo dnf install gcc gcc-c++ gcc-gfortran scons cmake
#   sudo dnf install boost-devel boost-python3 boost-python3-devel
#   sudo dnf install hdf5-devel suitesparse suitesparse-devel lapack-devel zlib-devel metis-devel

from scons.templates.rhel_options import *

# Disable MPI
mpi = 'none'
mpi4py = False
