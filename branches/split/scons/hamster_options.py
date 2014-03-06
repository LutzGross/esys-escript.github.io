
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

escript_opts_version = 201
#prefix = '/usr/local'
#build_dir = 'build'
#cxx = 'g++-4.8'
#cc_flags = ''
cc_optim = '-O3 -march=native'
#cc_debug = '-g'
cxx_extra = '-Wextra -Wno-unused-parameter -g'
#ld_extra = ''
#werror = False
#debug = True
verbose = True
openmp = True
#omp_flags = '-fopenmp'
#omp_ldflags = '-fopenmp'
mpi = 'OPENMPI'
mpi_prefix = '/usr/lib/openmpi'
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
#boost_prefix = '/usr/local'
boost_libs = ['libboost_python-py27']
#cppunit_prefix = '/usr/local'
#cppunit_libs = ['cppunit']
netcdf = True
#netcdf_prefix = '/usr/local'
#netcdf_libs = ['netcdf_c++', 'netcdf']
parmetis = True
#parmetis_prefix = '/usr/local'
#parmetis_libs = ['parmetis', 'metis']
#papi = True
#papi_prefix = '/usr/local'
#papi_libs = ['papi']
#papi_instrument_solver = True
#mkl = True
#mkl_prefix = '/usr'
#mkl_libs = ['mkl_solver', 'mkl_em64t', 'mkl_core', 'guide', 'pthread']
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
#boomeramg = True
#boomeramg_prefix = '/usr/local'
#boomeramg_libs = ['HYPRE']
lapack = 'clapack'
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
lapack_libs = ['lapack_atlas']
silo = True
#silo_prefix = '/usr/local'
silo_libs = ['siloh5']
visit = False
visit_prefix = '/opt/visit/2.7.0b/linux-x86_64/libsim/V2'
#visit_libs = ['simV2']
#build_shared = True


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

#vsl_random = True
#sys_libs = []
#env_export = []
#tools_names = ['default']
#iknowwhatimdoing = False
#forcelazy = 'leave_alone'
#forcecollres = 'leave_alone'

