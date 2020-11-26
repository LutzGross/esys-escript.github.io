
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

build_debug=False
build_debug=True

openmp=False
openmp=True

if build_debug == True:
	boost_prefix='/usr/local/boost.1.74.0'
	boost_libs='boost_python38'
	build_dir='build_debug'
	cxx_extra="-O0 -g -pg -fdiagnostics-color -Wno-implicit-int-float-conversion"
	cxx_extra+=" -DOXLEY_ENABLE_DEBUG"
	debug=True
	ld_extra='-L/usr/lib/openmpi/'
	if openmp is True:
		cxx_extra+=" -fopenmp"
		trilinos_prefix='/usr/local/trilinos_nompi'
	else:
		trilinos_prefix='/usr/local/trilinos_noomp'
	paso=False
	pythoncmd="/usr/bin/python3"
	pythonlibname="python3.8"
	pythonlibpath="/usr/lib"
	pythonincpath="/usr/include/python3.8"
	umfpack=True
	umfpack_prefix=['/usr/include/','/usr/lib']
	trilinos=True
	verbose=True
	werror=False
else:
	boost_prefix='/usr/local/boost.1.74.0'
	boost_libs='boost_python38'
	build_dir='build_normal'
	cxx_extra="-O3 -funroll-loops -fdiagnostics-color"
	debug=False
	ld_extra='-L/usr/lib/openmpi/'
	if openmp is True:
		cxx_extra+=" -fopenmp"
		trilinos_prefix='/usr/local/trilinos_nompi'
	else:
		trilinos_prefix='/usr/local/trilinos_noomp'
	paso=False
	pythoncmd="/usr/bin/python3"
	pythonlibname="python3.8"
	pythonlibpath="/usr/lib"
	pythonincpath="/usr/include/python3.8"
	umfpack=True
	umfpack_prefix=['/usr/include/','/usr/lib']
	trilinos=True
	verbose=True
	werror=False


