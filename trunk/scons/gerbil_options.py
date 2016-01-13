
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from templates.wheezy_options import *

cc_optim = '-O3 -march=native'
cxx_extra = '-Wextra -Wno-unused-parameter -g'
cuda = True
nvccflags = "-ccbin=g++ -arch=sm_30 -DBOOST_NOINLINE='__attribute__((noinline))'"
debug = False
#debug = True
if debug:
    cxx_extra += ' -fno-omit-frame-pointer -fsanitize=address'
    ld_extra = '-fsanitize=address'
verbose = True
boost_libs = ['libboost_python-py27']
mpi = 'OPENMPI'
parmetis = True
umfpack = True
lapack = 'clapack'
silo = True

#domains = 'ripley,speckley'
#longindices = True

