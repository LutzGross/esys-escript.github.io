# #############################################################################
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

compressed_files = True

# python3-gdal no longer supports 'import gdal'. Esys needs to be updated
gdal = True

escript_opts_version = 203
cxx_extra = ''          #To allow string addition later
#cxx_extra = '-Wno-literal-suffix'
openmp = True
#mpi = 'OPENMPI'

# sympy currently not supported as Sympy > 1.2 not supported

cxx_extra +=' -Wno-stringop-truncation'
werror = False

import os
import platform

mpil='unknown'
with open('/usr/share/mpi-default-dev/debian_defaults') as f:
  for line in f.readlines():
    bits = line.split('=')
    if bits[0] == "ARCH_DEFAULT_MPI_IMPL":
      mpil = bits[1].strip()
mpi=mpil.upper()

from subprocess import check_output
arch = check_output(['dpkg-architecture','-qDEB_HOST_MULTIARCH']).strip().decode()
libdir = '/usr/lib/' + arch
incdir = '/usr/include/' + arch

d_mpi_path = libdir + '/' + mpil
mpi_prefix = os.path.split(os.path.realpath(d_mpi_path))[0]
mpi_prefix = d_mpi_path
mpi_libs = ['mpi']

netcdf = 4
netcdf_prefix = [ '/usr/include', libdir]

# CppUnit library/libraries to link against
cppunit_prefix = [ '/usr/include', libdir]
cppunit_libs = ['cppunit']

umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', libdir]
umfpack_libs = ['umfpack', 'blas', 'amd']

lapack_prefix = [incdir , libdir ]

silo = True
silo_libs = ['siloh5']
dudley_assemble_flags = '-funroll-loops'

# Makefile.export.Tpetra currently missing on Debian
#if platform.machine() in ('x86_64','aarch64', 'ppc64el', 'ppc64', 'ppc64le', 's390x'):
#    trilinos = True
#    trilinos_prefix=[ '/usr/include/trilinos', libdir]
trilinos = False

pythoncmd = '/usr/bin/python3'

import sys

# this can be used by options files importing us
boost_py3_prefix = [ '/usr/include' , libdir ]
boost_py3_libs = [f'boost_python{sys.version_info.major}{sys.version_info.minor}']

boost_libs = boost_py3_libs
boost_prefix = boost_py3_prefix

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

if platform.architecture()[0] == '32bit':
    cxx_extra += ' -Wno-maybe-uninitialized '

if platform.machine() in ('ppc', 'ppc64', 'ppc64le', 's390x'):
    cxx_extra += ' -Wno-strict-overflow -Wno-error=strict-overflow '

mathjax_path='file:///usr/share/javascript/mathjax/MathJax.js'
