
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

openmp = True
umfpack = True
silo = True
# mpi = 'OPENMPI'
# verbose = True
# debug = True
# trilinos = True
# paso = False
# parmetis = True
visit = True
#werror = False
# cxx = 'clang++'
# cxx='/usr/bin/clang++'

python = 3

##############################################################################
escript_opts_version = 203

import os
import subprocess

cxx_extra += " -fmessage-length=80 -fdiagnostics-color=always "
if trilinos is True:
  cxx_extra += "  -Wno-deprecated-declarations -Wno-unused-variable "

netcdf = 4
mpi_libs = ['mpi_cxx', 'mpi']
parmetis_libs = ['parmetis', 'metis']
silo_libs = ['siloh5', 'hdf5_cpp']
umfpack_libs = ['umfpack', 'blas', 'amd']

lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
d_mpi_path = '/usr/include/openmpi'
mpi_prefix = os.path.split(os.path.realpath(d_mpi_path))[0]
parmetis_prefix = ['/usr/include','/usr/lib']
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']

visit_prefix = ['/usr/local/visit/2.13.2/linux-x86_64/libsim/V2/include/','/usr/local/visit/2.13.2/linux-x86_64/libsim/V2/lib/']

if cxx=='clang++':
  trilinos_prefix =['/usr/local/trilinos_clang/include/','/usr/local/trilinos_clang/lib/']
elif mpi=='OPENMPI':
  trilinos_prefix =['/usr/local/trilinos_mpi/include/','/usr/local/trilinos_mpi/lib/']
else:
  trilinos_prefix =['/usr/local/trilinos/include/','/usr/local/trilinos/lib/']

if trilinos == True:
  print("Meerkat config: Using trilinos libraries %s" % trilinos_prefix) 

 
boost_prefix='/usr/local'
p = subprocess.Popen(["ld","--verbose"], stdout=subprocess.PIPE)
out,err = p.communicate()
spath = [x[13:-3] for x in out.split() if 'SEARCH_DIR' in x]
p2name = ''
p3name = ''
for name in spath:
  try:
    l=os.listdir(name)
    p2res=[x for x in l if x.startswith('libboost_python2') and x.endswith('.so')]
    p3res=[x for x in l if x.startswith('libboost_python3') and x.endswith('.so')]
    if len(p2name)==0 and len(p2res)>0:
      p2name=p2res[-1]
    if len(p3name)==0 and len(p3res)>0:
      p3name=p3res[-1]
  except OSError:
    pass

# boost-python library/libraries to link against
if python == 2:
  boost_libs = [p2name[3:-3]]
  pythoncmd = '/usr/bin/python'
else:
  boost_libs = [p3name[3:-3]]
  pythoncmd = '/usr/bin/python3'

print("Meerkat config: Linking with %s" % boost_libs) 

#boost_libs = [p2name[3:-3], 'boost_numpy27']
# boost_libs = ['boost_python27', 'boost_numpy27']

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
