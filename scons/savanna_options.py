
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

# module load intel intel-mpi silo/4.10.2 netcdf/4.6.1 hdf5/1.8.14 gdal/2.1.0 boost/1.61.0 boost/1.73.0 cppunit/1.13.2 scons/2.5.0 python/3.6.9 visit/2.10.0 gmsh/3.0.6-intel

escript_opts_version = 203
import os
MKLROOT = os.environ['MKLROOT']
I_MPI_ROOT = os.environ['I_MPI_ROOT']

boost_prefix = '/sw/libs/boost/1.73.0'
boost_libs = ['boost_python36','boost_numpy36','boost_iostreams']
build_trilinos='never'
build_dir='build.new'
cppunit_prefix = '/sw/libs/cppunit/1.13.2'
cxx='mpiicpc'
cxx_extra=' -I'+os.path.join(MKLROOT, 'include ')
# cxx_extra+=' -march=native -lto '
domains=['finley','oxley','ripley','speckley']
gmsh=1
lapack=True
lapack_prefix=[os.path.join(MKLROOT,'include'),os.path.join(MKLROOT, 'lib','intel64')]
launcher = "srun --nodes=%n --ntasks=%N --ntasks-per-node=%p --cpus-per-task=%t --cpu_bind=quiet %b"
ld_extra = '-ipo-separate -shared-intel -L'+os.path.join(MKLROOT, 'lib','intel64')+' -L/sw/libs/hdf5/1.8.14/lib'
ld_extra += ' -wd11021 -wd3180 -wd161 -Wno-unknown-pragmas -wd11012 ' 
ld_extra += ' -L/sw/libs/hdf5/1.8.14/lib/ -lhdf5 -L'+os.path.join(MKLROOT, 'lib','intel64')
mkl = True
_mklroot=MKLROOT or '/sw/intel/compilers_and_libraries_2018.1.163/linux/mkl'
mkl_prefix = ['%s/include'%_mklroot, '%s/lib/intel64'%_mklroot]
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core']
mpi='INTELMPI'
mpi_prefix = I_MPI_ROOT or '/sw/intel/impi/2018.1.163'
mpi_prefix += '/intel64'
netcdf=4
netcdf_prefix='/sw/libs/netcdf/4.6.1'
netcdf_libs = ['netcdf', 'netcdf_c++4']
openmp=1
parmetis = True
parmetis_prefix = '/sw/libs/parmetis/4.0.3-impi'
parmetis_libs = ['parmetis']
paso=0
postlaunch = ""
prelaunch = ""
pythoncmd="/sw/apps/python/3.6.9/bin/python3"
pythonlibname='python3'
silo = True
silo_prefix = '/sw/libs/silo/4.10.2'
silo_libs=['siloh5']
build_trilinos='make'
trilinos_make=os.getcwd()+'/scons/savanna_trilinosmake.sh'
print(os.getcwd())
use_p4est=0
umfpack=0
use_p4est=1
verbose=0
visit=1
visit_prefix='/sw/apps/visit/2.10.0/linux-x86_64/libsim/V2/'
werror=0

