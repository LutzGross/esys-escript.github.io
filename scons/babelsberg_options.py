
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from templates.ubuntu_options import *

# MANDATORY: Explicitly set version to ensure SCons recognizes it
escript_opts_version = 203

mpi = "none"

# LAPACK configuration - uses LAPACKE (modern C interface)
lapack = 'auto'  # Auto-detect LAPACKE
lapack_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
lapack_libs = ['lapacke']
