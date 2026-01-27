
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

# No-MPI configuration for Linux Mint.
# Imports the full mint_options.py and disables MPI.

from scons.templates.mint_options import *

# Disable MPI
mpi = 'none'
mpi4py = False

# ParMETIS requires MPI
parmetis = False
