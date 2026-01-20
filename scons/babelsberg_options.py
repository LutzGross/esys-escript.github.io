
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
build_trilinos = 'make'  # Use smart dependency checking

# Include oxley with sc_MPI_COMM_WORLD for p4est/p8est integration
domains = ['finley', 'ripley', 'speckley', 'oxley']

# Override cxx_extra to avoid duplicate optimization flags with MPI
cxx_extra = ['-fdiagnostics-color=always', '-fstack-protector-strong', '-Wformat', '-Werror=format-security']
