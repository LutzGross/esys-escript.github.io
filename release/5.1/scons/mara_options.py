
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

from templates.sid_py3_mpi_options import *

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
parmetis = True
trilinos = True
trilinos_prefix = '/opt/trilinos_hybrid_eti'
umfpack = True
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']
silo = True
visit = False
visit_prefix = '/opt/visit/2.7.0b/linux-x86_64/libsim/V2'
#visit_libs = ['simV2']
launcher = "mpirun ${AGENTOVERRIDE} ${EE} --map-by node:pe=%t -bind-to none -np %N %b"

#longindices = True
#cxx_extra += ' -Wconversion'
#lapack = 'none'
#parmetis = False

