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
from templates.jessie_mpi_options import *
cc_optim = '-O3 -march=native'
verbose = True
parmetis = True
umfpack = True
silo = True
trilinos=True
trilinos_prefix="/opt/trilinos_hybrid"
werror=False


parmetis_prefix='/usr/local'
#longindices=True
