
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

build_debug=False
#build_debug=True

openmp=False
openmp=True

use_paso=False

if build_debug == True:
	#boost_prefix='/usr/local/boost.1.74.0'
	boost_libs='boost_python39'
	build_dir='build_debug'
	cxx_extra="-O0 -g -pg -fdiagnostics-color -Wno-implicit-int-float-conversion"
	cxx_extra+=" -Wno-deprecated-declarations"
	# cxx_extra+=" -DOXLEY_ENABLE_DEBUG"
	# cxx_extra+=" -DOXLEY_ENABLE_DEBUG_NODES"
	debug=True
	ld_extra='-L/usr/lib/openmpi/'
	if openmp is True:
		cxx_extra+=" -fopenmp"
		trilinos_prefix='/usr/local/trilinos_nompi'
	else:
		trilinos_prefix='/usr/local/trilinos_noomp'
	paso=use_paso
	pythoncmd="/usr/bin/python3"
	pythonlibname="python3.9"
	pythonlibpath="/usr/lib/x86_64-linux-gnu/"
	pythonincpath="/usr/include/python3.9"
	umfpack=True
	umfpack_prefix=['/usr/include/suitesparse','/usr/lib/x86_64-linux-gnu']
	trilinos=True
	verbose=True
	werror=False
else:
	#boost_prefix='/usr/local/boost.1.74.0'
	boost_libs='boost_python39'
	build_dir='build_normal'
	cxx_extra="-O2 -funroll-loops -march=native -fdiagnostics-color"
	debug=False
	ld_extra='-L/usr/lib/openmpi/'
	if openmp is True:
		cxx_extra+=" -fopenmp"
		trilinos_prefix='/usr/local/trilinos_nompi'
	else:
		trilinos_prefix='/usr/local/trilinos_noomp'
	paso=use_paso
	pythoncmd="/usr/bin/python3"
	pythonlibname="python3.9"
	pythonlibpath="/usr/lib/x86_64-linux-gnu/"
	pythonincpath="/usr/include/python3.9"
	umfpack=True
	umfpack_prefix=['/usr/include/suitesparse','/usr/lib/x86_64-linux-gnu']
	trilinos=True
	verbose=True
	werror=False

cxx_extra+=" -Wno-maybe-unitialized"
