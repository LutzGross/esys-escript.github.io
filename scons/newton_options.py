
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
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

from templates.trusty_options import *

mpi = 'OPENMPI'
boost_libs = ['boost_python-py27']
parmetis = True
umfpack = True
lapack = 'clapack'
lapack_libs = ['lapack']
silo = True
silo_libs = ['siloh5', 'hdf5']
launcher = "mpirun --gmca mpi_warn_on_fork 0 ${EE} -np %N %b"

