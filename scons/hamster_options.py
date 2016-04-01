
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

from templates.wheezy_options import *
cc_optim = '-O3 -march=native'
#cc_debug = "-g3 -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK -D_GLIBCXX_DEBUG -fno-omit-frame-pointer -fsanitize=address --param=max-vartrack-size=90000000"
cxx_extra = '-Wextra -Wno-unused-parameter'
verbose = False
mpi = 'OPENMPI'
#ld_extra = '-fsanitize=address'
boost_libs = ['boost_python-py27']
parmetis = True
umfpack = True
lapack = 'clapack'
silo = True

