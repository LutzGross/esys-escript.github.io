
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# This is a template configuration file for escript on Fedora Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
pythoncmd='/usr/bin/python3'
pythonlibpath = ['/usr/lib64']
pythonlibname = ['python3.10']
pythonincpath = ['/usr/include/python3.10']
boost_libs = ['boost_python310']
boost_prefix=['/usr/include','/usr/lib64']
disable_boost_numpy=False
umfpack=True
umfpack_prefix=['/usr/include/','/usr/lib64']
werror=0
