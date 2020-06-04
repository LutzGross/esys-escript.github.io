##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

#from templates.jessie_mpi_options import *
#from templates.jessie_py3_options.py import *
from templates.buster_py3_options import *

verbose = True
silo = True
#boost_libs = ['boost_python']
umfpack = True
werror=False
#pythoncmd = '/home/gross/anaconda3/bin/python3'
pythoncmd="/usr/bin/python3"
pythonlibname = 'python3.6m'
#pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = '/usr/include/python3.6m'

silo_libs = ['siloh5', 'hdf5']
silo_prefix=[ '/usr/include' , '/usr/lib/x86_64-linux-gnu/hdf5/serial' ]
#netcdf_prefix=[ '/usr/local/2.13.0/linux-x86_64/include/netcdf/include/' ]


netcdf = 3



trilinos=True
trilinos_prefix="/opt/trilinos"

openmp = True 
parmetis = False

mpi = 'no' # 'OPENMPI'

#pythoncmd = 'python3'
#boost_libs = boost_py3_libs

