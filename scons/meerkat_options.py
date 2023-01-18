
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

debug=0
debug=1

openmp=False
openmp=True

paso=1
build_trilinos=0
trilinos=0
trilinos_prefix='~/Documents/escript6.0/escript-trilinos'

if debug == True:
	boost_libs='boost_python39'
	build_dir='build_debug'
	cxx_extra="-O0 -g -pg -fdiagnostics-color -Wno-implicit-int-float-conversion"
	debug=True
	ld_extra='-L/usr/lib/openmpi/'
	pythoncmd="/usr/bin/python3"
	pythonlibname="python3.9"
	pythonlibpath="/usr/lib/x86_64-linux-gnu/"
	pythonincpath="/usr/include/python3.9"
	umfpack=True
	umfpack_prefix=['/usr/include/suitesparse','/usr/lib/x86_64-linux-gnu']
	verbose=True
	werror=False
else:
	boost_libs='boost_python39'
	build_dir='build_normal'
	cxx_extra="-O3 -funroll-loops -fdiagnostics-color"
	ld_extra='-L/usr/lib/openmpi/'
	pythoncmd="/usr/bin/python3"
	pythonlibname="python3.9"
	pythonlibpath="/usr/lib/x86_64-linux-gnu/"
	pythonincpath="/usr/include/python3.9"
	umfpack=True
	umfpack_prefix=['/usr/include/suitesparse','/usr/lib/x86_64-linux-gnu']
	verbose=True
	werror=False

cxx_extra+=" -Wno-maybe-unitialized"

# use_sympy=False