

##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


# ADD TO bash_profile

## module load libraries/boost/1.67-gnu-7.3
## module load softwares/python/3.6-anaconda-5.0.1
## module load compilers/intel/2018.2.199
## module load libraries/mkl/intel-2018.2.199
## #module load libraries/parmetis/4.0.3-intel-2018.2.199
## module load compilers/gnu/7.3
## eval $loadIntelLibs
## export ESYSROOT="/home/lgross/escript"
## export PATH=$ESYSROOT/bin:${PATH}
## export LD_LIBRARY_PATH=${ESYSROOT}/lib:${LD_LIBRARY_PATH}
## export PYTHONPATH=${ESYSROOT}:${PYTHONPATH}



# this is for sc.npad.imd.ufrn.br 

MKLROOT=None
try:
    import os
    #if not 'escript/dev-deps' in os.environ['LOADEDMODULES'].split(':'):
    #    print("WARNING: The escript/dev-deps module does not appear to be loaded!")
    MKLROOT = os.environ['MKLROOT']
    I_MPI_ROOT = os.environ['I_MPI_ROOT']
except:
    MKLROOT = None
    I_MPI_ROOT = None

escript_opts_version = 203

#cxx_extra = '-sox -I/opt/npad/shared/softwares/python/3.6.3-gnu-4.8_sh_lib/lib/python3.6/site-packages/numpy/core/include/'
cxx_extra = '-sox -I/opt/npad/shared/softwares/python/3.6-anaconda-5.0.1/pkgs/numpy-1.13.3-py36ha12f23b_0/lib/python3.6/site-packages/numpy/core/include'

#ld_extra = '-ipo-separate -shared-intel -L/sw/libs/hdf5/1.8.14/lib'
#ld_extra = '-shared-intel -L/sw/libs/hdf5/1.8.14/lib'
ld_extra = '-shared-intel'
ld_extra += ' -wd11021 '  #silence icpc warnings about symbols ipo can't see

werror = False
verbose = True
openmp = True

pythoncmd = "python3"
pythonlibname ="libpython3.6m.a"
pythonlibpath="/opt/npad/shared/softwares/python/3.6-anaconda-5.0.1/lib"

# mpi = 'NONE'
mpi = 'INTELMPI'
mpi_prefix = I_MPI_ROOT or '/opt/intel/parallel_studio_xe_2018/compilers_and_libraries_2018.2.199/linux/mpi'
mpi_prefix += '/intel64'
#cuda = True
#nvccflags = "-arch=sm_35 -ccbin=icpc -DBOOST_NOINLINE='__attribute__((noinline))'"
#cuda_prefix = ['/sw/libs/cuda/7.5/include', '/sw/libs/cuda/7.5/lib64']

boost_prefix = '/opt/npad/shared/libraries/boost/1.67-gnu-7.3'
boost_libs = ['boost_python36']

#cppunit_prefix = '/sw/libs/cppunit/1.13.2'
netcdf = False
netcdf_prefix = '/opt/npad/shared/libraries/netcdf-c/4.3.0-intel-2018.2.199'
netcdf_libs = ['netcdf_c++', 'netcdf', 'hdf5']

parmetis = False
parmetis_prefix = '/opt/npad/shared/libraries/parmetis/4.0.3-intel-2018.2.199/'
parmetis_libs = ['parmetis']

trilinos = False
trilinos_prefix = '/sw/libs/trilinos/snapshot-hybrid-eti'

mkl = True
_mklroot=MKLROOT
mkl_prefix = ['%s/include'%_mklroot, '%s/lib/intel64'%_mklroot]
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']
#boomeramg = True
#boomeramg_prefix = '/sw/libs/hypre/2.0.0'
#boomeramg_libs = ['HYPRE']
#boomeramg_libs = ['HYPRE_IJ_mv', 'HYPRE_krylov', 'HYPRE_parcsr_ls']
silo = False
silo_prefix = '/sw/libs/silo/4.10.2'
silo_libs = ['siloh5', 'hdf5']
visit = False
# visit_prefix = '/sw/apps/visit/2.7.0/linux-x86_64/libsim/V2'

#prelaunch = ""
launcher = "srun --nodes=%n --ntasks=%N --ntasks-per-node=%p --cpus-per-task=%t --cpu_bind=quiet %b"
#postlaunch = ""

env_export = ['INTEL_LICENSE_FILE']

tools_names = [('intelc',{'topdir':'/sw/intel/compilers_and_libraries_2018.1.163'})]

# uncomment the following four options to build with mpt (check modules!)
#build_dir = 'buildmpt'
#mpi = 'MPT'
#mpi_prefix = '/opt/sgi/mpt/mpt-2.10'

