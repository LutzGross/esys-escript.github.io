
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
escript_opts_version = 2032


# Enable MPI - auto-detect flavour from mpi4py
mpi = 'auto'

# Enable mpi4py support
mpi4py = True

# Enable parmetis
parmetis = True
parmetis_prefix = ['/usr/include/parmetis', '/usr/lib/x86_64-linux-gnu']
parmetis_libs = ['parmetis']

# Build Trilinos with MPI support
build_trilinos = 'always'

# Exclude oxley since it requires Trilinos
domains = ['finley', 'ripley', 'speckley']

# Override cxx_extra to avoid duplicate optimization flags with MPI
cxx_extra = ['-fdiagnostics-color=always', '-fstack-protector-strong', '-Wformat', '-Werror=format-security']
