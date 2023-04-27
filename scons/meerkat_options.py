
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

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203

#cxx='/usr/bin/clang++'

debug=False
# debug=True

# openmp=False
openmp=True

boost_libs='boost_python39'
build_dir='build_normal'
cxx_extra=" -fdiagnostics-color"

# cxx_extra+=" -p -g -pg -ggdb -O0"

domains=["finley"]

ld_extra='-L/usr/lib/openmpi/'
paso=False
pythoncmd="/usr/bin/python3"
pythonlibname="python3.9"
pythonlibpath="/usr/lib/x86_64-linux-gnu/"
pythonincpath="/usr/include/python3.9"

cxx_extra+=" -Wno-maybe-unitialized"

# paso=1
# build_trilinos=0
build_trilinos='never'
trilinos=1
trilinos_prefix="/usr/local/trilinos.mpi.13"
mpi='OPENMPI'
mpi_prefix='/usr/lib/x86_64-linux-gnu/openmpi/'

verbose=1
werror=0



