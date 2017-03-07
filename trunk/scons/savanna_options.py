
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

MKLROOT=None
try:
    import os
    if not 'escript/dev-deps' in os.environ['LOADEDMODULES'].split(':'):
        print("WARNING: The escript/dev-deps module does not appear to be loaded!")
    MKLROOT = os.environ['MKLROOT']
    I_MPI_ROOT = os.environ['I_MPI_ROOT']
except:
    pass

escript_opts_version = 203

cxx_extra = '-sox -I/sw/pymodules/3.5/scipy-0.17.1-haswell/lib/python3.5/site-packages/numpy-1.11.0-py3.5-linux-x86_64.egg/numpy/core/include'

ld_extra = '-ipo-separate -shared-intel -L/sw/libs/hdf5/1.8.14/lib'
ld_extra += ' -wd11021 '  #silence icpc warnings about symbols ipo can't see

werror = False
verbose = True
openmp = True

pythoncmd = "python3"

mpi = 'INTELMPI'
mpi_prefix = I_MPI_ROOT or '/sw/intel/impi/2017.0.098'
mpi_prefix += '/intel64'
#cuda = True
nvccflags = "-arch=sm_35 -ccbin=icpc -DBOOST_NOINLINE='__attribute__((noinline))'"
cuda_prefix = ['/sw/libs/cuda/7.5/include', '/sw/libs/cuda/7.5/lib64']

boost_prefix = '/sw/libs/boost/1.61.0'
boost_libs = ['boost_python3']
cppunit_prefix = '/sw/libs/cppunit/1.13.2'
netcdf = True
netcdf_prefix = '/sw/libs/netcdf/4.1.3'
netcdf_libs = ['netcdf_c++', 'netcdf', 'hdf5']
parmetis = True
parmetis_prefix = '/sw/libs/parmetis/4.0.3-impi'
parmetis_libs = ['parmetis']
trilinos = True
trilinos_prefix = '/sw/libs/trilinos/snapshot-hybrid-eti'

mkl = True
_mklroot=MKLROOT or '/sw/intel/compilers_and_libraries_2017.0.098/linux/mkl'
mkl_prefix = ['%s/include'%_mklroot, '%s/lib/intel64'%_mklroot]
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']
#boomeramg = True
boomeramg_prefix = '/sw/libs/hypre/2.0.0'
boomeramg_libs = ['HYPRE']
#boomeramg_libs = ['HYPRE_IJ_mv', 'HYPRE_krylov', 'HYPRE_parcsr_ls']
silo = True
silo_prefix = '/sw/libs/silo/4.10.2'
silo_libs = ['siloh5', 'hdf5']
visit = False
visit_prefix = '/sw/apps/visit/2.7.0/linux-x86_64/libsim/V2'

prelaunch = ""
launcher = "srun --nodes=%n --ntasks=%N --ntasks-per-node=%p --cpus-per-task=%t --cpu_bind=quiet %b"
postlaunch = ""

env_export = ['INTEL_LICENSE_FILE']

tools_names = [('intelc',{'topdir':'/sw/intel/compilers_and_libraries_2017.0.098'})]

# uncomment the following four options to build with mpt (check modules!)
#build_dir = 'buildmpt'
#mpi = 'MPT'
#mpi_prefix = '/opt/sgi/mpt/mpt-2.10'
#parmetis_prefix = '/sw/libs/parmetis/4.0.3-mpt'

