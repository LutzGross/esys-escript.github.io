
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

mpi_no_host=True

cxx_extra = '-Wno-deprecated-declarations'
openmp = True
boost_libs = ['boost_python-py27']
mpi = 'OPENMPI' 
mpi_prefix = '/usr/lib/x86_64-linux-gnu/openmpi'
#mpi_prefix = '/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/'
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
netcdf = True
#umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
silo = True
silo_libs = ['siloh5']
dudley_assemble_flags = '-funroll-loops'

from site_init import getdebbuildflags
# Now we add the debian build flags
debstuff = getdebbuildflags()
if len(debstuff) > 0:
  print("Building with the following additional flags from debian: "+str(debstuff))
for i in debstuff:
  k=i[0]
  v=i[1]
  try:
    exec(k+"+=' "+v+"'")
  except NameError:   
    exec(k+"='"+v+"'")

trilinos=True

trilinos_prefix='/opt/trilinos_hybrid'
parmetis=True

parmetis_prefix='/usr/local'
