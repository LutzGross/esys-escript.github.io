
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

from templates.stretch_py3_mpi_options import *

#copied this from guineapig originally
# disabled until the boost issue is fixed.
#cuda = True
#nvccflags = "-ccbin=g++-4.9 -arch=sm_30 -DBOOST_NOINLINE='__attribute__((noinline))'"

parmetis = True
umfpack = True
silo = True
trilinos = True
trilinos_prefix = '/usr/local/trilinos/post12.12-1'
#cxx_extra += " -Wextra -Wno-deprecated-declarations -Wno-unused-parameter"
launcher = "mpirun ${AGENTOVERRIDE} ${EE} --mca io romio314 --oversubscribe --map-by node:pe=%t -bind-to none -np %N %b"

