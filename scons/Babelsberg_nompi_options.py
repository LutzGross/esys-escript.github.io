
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

# No-MPI test configuration for release testing
from templates.ubuntu_options import *

# MANDATORY: Explicitly set version to ensure SCons recognizes it
escript_opts_version = 2032

# Disable MPI
mpi = 'none'
mpi4py = False

# Disable Trilinos (requires MPI for full functionality)
trilinos = True
build_trilinos = 'make'
#build_trilinos = 'always'

# Domains without oxley (matching release tarball)
domains = ['finley', 'ripley', 'speckley']

# Override cxx_extra
cxx_extra = ['-fdiagnostics-color=always', '-fstack-protector-strong', '-Wformat', '-Werror=format-security']
