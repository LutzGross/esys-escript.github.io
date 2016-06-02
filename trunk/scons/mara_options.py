
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

escript_opts_version = 203
#cuda = True
cc_optim = '-O3 -march=native'
cc_debug = "-g3 -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK -D_GLIBCXX_DEBUG -fno-omit-frame-pointer" #-fsanitize=address 
cxx_extra = '-Wextra -Wno-unused-parameter -Wno-deprecated-declarations -g -fdiagnostics-color'
nvccflags = "-arch=sm_30 -DBOOST_NOINLINE='__attribute__((noinline))'"
#ld_extra = ''
#werror = False
#debug = True
#ld_extra = '-fsanitize=address'
verbose = True
openmp = True
mpi = 'OPENMPI'
mpi_prefix = '/usr/lib/openmpi'
mpi_libs = ['mpi_cxx', 'mpi']
boost_libs = ['libboost_python-py27']
netcdf = True
parmetis = True
trilinos = True
trilinos_prefix = '/opt/trilinos'
#papi = True
#papi_prefix = '/usr/local'
#papi_libs = ['papi']
#papi_instrument_solver = True
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
silo = True
silo_libs = ['siloh5']
visit = False
visit_prefix = '/opt/visit/2.7.0b/linux-x86_64/libsim/V2'
#visit_libs = ['simV2']

#longindices = True
#cxx_extra += ' -Wconversion'
#lapack = 'none'
#parmetis = False

