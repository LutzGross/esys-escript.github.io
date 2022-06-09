
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

clang=False
# clang=True

build_debug=False
build_debug=True

openmp=False
openmp=True

trilinos_debug=False
# trilinos_debug=True

use_paso=False
# use_paso=True

verbose=0
# verbose=1
# verbose=2

##############################################################################

escript_opts_version=203
boost_libs='boost_python39'	
cxx_extra=""
paso=use_paso
pythoncmd="/usr/bin/python3"
pythonlibname="python3.9"
pythonlibpath="/usr/lib/x86_64-linux-gnu/"
pythonincpath="/usr/include/python3.9"
silo=True
silo_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu']
silo_libs=['siloh5','hdf5_serial']
umfpack=True
umfpack_prefix=['/usr/include/suitesparse','/usr/lib/x86_64-linux-gnu']
trilinos=True

if openmp is True:
	cxx_extra+=" -fopenmp"
	trilinos_prefix='/usr/local/trilinos_nompi'
else:
	trilinos_prefix='/usr/local/trilinos_noomp'

if clang is True:
	cxx='/usr/bin/clang++'
	cxx_extra+=" -Wno-implicit-const-int-float-conversion -Wno-overloaded-virtual"
	cxx_extra+=" -Wno-unknown-warning-option"
	cxx_extra+=" -Wno-non-c-typedef-for-linkage"
	cxx_extra+=" -Wno-dtor-name"
	cxx_extra+=" -Wno-overloaded-virtual"
	if build_debug is True:
		build_dir='build_clang'
	else:
		build_dir='build_clang_norm'
else:
	cxx_extra+=" -Wno-deprecated-declarations"
	if build_debug is True:
		build_dir='build_debug'
	else:
		build_dir='build_normal'

if build_debug is True:
	debug=True
	cxx_extra+=" -O0 -g -pg -fdiagnostics-color -Wno-implicit-int-float-conversion"
	cxx_extra+=" -Wno-unused-function"
	cxx_extra+=" -Wno-unused-variable"
	cxx_extra+=" -Wno-unused-but-set-variable"
	build_full=1
	werror=False
else:
	debug=False
	cxx_extra+=" -O2 -march=native -fdiagnostics-color"
	cxx_extra+="-funroll-loops -flto "

if trilinos_debug is True:
	trilinos_prefix='/usr/local/trilinos_debug'

if verbose==2:
	verbose=1
	cxx_extra+=" --verbose"

cxx_extra+=" -Wno-maybe-unitialized"
