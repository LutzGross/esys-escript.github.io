
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

# No-MPI configuration for Debian/Ubuntu Linux.
# Imports the full debian_options.py and disables MPI.
#
# Prerequisites (MPI packages not required):
#   sudo apt-get install python3-dev python3-numpy python3-scipy python3-sympy
#   sudo apt-get install scons cmake g++ gfortran
#   sudo apt-get install libboost-python-dev libboost-numpy-dev libboost-random-dev libboost-iostreams-dev
#   sudo apt-get install libnetcdf-dev libsilo-dev libhdf5-serial-dev
#   sudo apt-get install libsuitesparse-dev liblapacke-dev libmumps-seq-dev
#   sudo apt-get install libmetis-dev zlib1g-dev

from debian_options import *

# Disable MPI
mpi = 'none'
mpi4py = False

# ParMETIS requires MPI
parmetis = False

# Don't build Trilinos (optional, but faster build without MPI)
build_trilinos = 'never'
trilinos = False
