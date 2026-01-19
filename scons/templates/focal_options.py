
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file in the repository root for contributors and development history
# https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS
#
##############################################################################

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203

openmp = True

umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

pythoncmd="/usr/bin/python3"

import subprocess
p=subprocess.run([pythoncmd,'-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1]
pythonlibname = 'python3.%s'%subversion
pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = '/usr/include/python3.%s'%subversion

boost_libs=['boost_python3%s'%subversion,'boost_numpy3%s'%subversion,'boost_random']
boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
