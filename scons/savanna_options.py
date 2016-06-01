
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

try:
    import os
    if not 'escript/dev-deps' in os.environ['LOADEDMODULES'].split(':'):
        print("WARNING: The escript/dev-deps module does not appear to be loaded!")
except:
    pass

escript_opts_version = 203

cxx_extra = '-sox -I/sw/pymodules/2.7/scipy-0.15.1-haswell/lib/python2.7/site-packages/numpy/core/include'

ld_extra = '-ipo-separate -shared-intel -L/sw/libs/hdf5/1.8.14/lib'
ld_extra += ' -wd11021 '  #silence icpc warnings about symbols ipo can't see

werror = False
verbose = True
openmp = True

mpi = 'INTELMPI'
mpi_prefix = '/sw/intel/impi/5.1.1.109/intel64'
#cuda = True
nvccflags = "-arch=sm_35 -ccbin=icpc -DBOOST_NOINLINE='__attribute__((noinline))'"
cuda_prefix = ['/sw/libs/cuda/7.5/include', '/sw/libs/cuda/7.5/lib64']

boost_prefix = '/sw/libs/boost/1.57.0'
boost_libs = ['boost_python']
cppunit_prefix = '/sw/libs/cppunit/1.13.2'
netcdf = True
netcdf_prefix = '/sw/libs/netcdf/4.1.3'
netcdf_libs = ['netcdf_c++', 'netcdf', 'hdf5']
parmetis = True
parmetis_prefix = '/sw/libs/parmetis/4.0.3-impi'
parmetis_libs = ['parmetis']

mkl = True
_mklroot='/sw/intel/compilers_and_libraries_2016.0.109/linux/mkl'
mkl_prefix = ['%s/include'%_mklroot, '%s/lib/intel64'%_mklroot]
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']
#boomeramg = True
boomeramg_prefix = '/sw/libs/hypre/2.0.0'
boomeramg_libs = ['HYPRE']
#boomeramg_libs = ['HYPRE_IJ_mv', 'HYPRE_krylov', 'HYPRE_parcsr_ls']
lapack = 'mkl'
lapack_prefix = mkl_prefix
lapack_libs = mkl_libs
silo = True
silo_prefix = '/sw/libs/silo/4.10.2'
silo_libs = ['siloh5', 'hdf5']
visit = False
visit_prefix = '/sw/apps/visit/2.7.0/linux-x86_64/libsim/V2'

prelaunch = ""
launcher = "srun --nodes=%n --ntasks=%N --ntasks-per-node=%p --cpus-per-task=%t --cpu_bind=quiet %b"
postlaunch = ""

env_export = ['INTEL_LICENSE_FILE']

tools_names = [('intelc',{'topdir':'/sw/intel/compilers_and_libraries_2016.0.109'})]

# uncomment the following four options to build with mpt (check modules!)
#build_dir = 'buildmpt'
#mpi = 'MPT'
#mpi_prefix = '/opt/sgi/mpt/mpt-2.10'
#parmetis_prefix = '/sw/libs/parmetis/4.0.3-mpt'

