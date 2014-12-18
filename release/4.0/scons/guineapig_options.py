
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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
cuda = True

nvccflags = "-ccbin=g++-4.8 -arch=sm_30 -DBOOST_NOINLINE='__attribute__((noinline))'"

mpi = 'OPENMPI'

boost_libs = ['boost_python-py27']

parmetis = True

umfpack = True

lapack = 'clapack'

silo = True

silo_libs = ['siloh5', 'hdf5_openmpi']

build_shared = False

launcher = "mpirun --gmca mpi_warn_on_fork 0 -x ${EE} -np %N %b"
