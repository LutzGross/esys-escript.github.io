
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
#cxx_extra = '-Wno-literal-suffix'
openmp = True
#mpi = 'OPENMPI'
umfpack = True

#pythoncmd='/usr/bin/python3.5'
pythoncmd='/usr/bin/python2.7'

import os

d_mpi_path = '/usr/include/openmpi'
mpi_prefix = os.path.split(os.path.realpath(d_mpi_path))[0]
mpi_libs = ['mpi_cxx', 'mpi']
netcdf = 4
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
#silo = True
silo_libs = ['siloh5', 'hdf5_openmpi']
dudley_assemble_flags = '-funroll-loops'


import subprocess
import os
p = subprocess.Popen(["ld","--verbose"], stdout=subprocess.PIPE)
out,err = p.communicate()
spath = [x[13:-3] for x in out.split() if 'SEARCH_DIR' in x]
p2name = ''
p3name = ''
for name in spath:
  try:
    l=os.listdir(name)
    p2res=[x for x in l if x.startswith('libboost_python-py2') and x.endswith('.so')]
    p3res=[x for x in l if x.startswith('libboost_python-py3') and x.endswith('.so')]
    if len(p2name)==0 and len(p2res)>0:
      p2name=p2res[-1]
    if len(p3name)==0 and len(p3res)>0:
      p3name=p3res[-1]
  except OSError:
    pass

# boost-python library/libraries to link against
boost_libs = [p2name[3:-3]]
#boost_libs = [p3name[3:-3]]

# this can be used by options files importing us
boost_py2_libs = [p2name[3:-3]]
boost_py3_libs = [p3name[3:-3]]

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

mathjax_path='/usr/share/javascript/mathjax/MathJax.js'
