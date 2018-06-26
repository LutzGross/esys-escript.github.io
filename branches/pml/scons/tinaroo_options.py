
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

# required modules in env:
#module load intel_xe/2016.2.181
#module load GCC/4.9.2
#module load intel_mpi/5.0.3.048

escript_opts_version = 203

prefix = '/30days/uqcalti1/prefix'

cc_optim = '-O3 -ftz -fno-alias -inline-level=2 -ipo -xCORE-AVX2 -funroll-loops'
cxx_extra = '-wd2259 -Wall -Wcheck -sox -I%s/lib/python3.5/site-packages/numpy/core/include'%prefix

#ld_extra = '-ipo-separate -shared-intel -qopt-report=5 -qopt-report-phase=vec,loop'
ld_extra = '-ipo-separate -shared-intel '
ld_extra += ' -wd11021 '  #silence icpc warnings about symbols ipo can't see

werror = False
verbose = True
openmp = True
omp_flags = '-qopenmp'
omp_ldflags = '-qopenmp'
mpi = 'INTELMPI'
mpi_prefix = '/sw/Licensed/intel/impi/5.0.3.048/intel64'
pythoncmd = 'python3.5'
boost_prefix = prefix
boost_libs = ['boost_python3']
netcdf = True
netcdf_prefix = prefix
netcdf_libs = ['netcdf_c++', 'netcdf', 'hdf5']
parmetis = False
mkl = True
_mklroot='/sw/Licensed/intel/compilers_and_libraries_2016.2.181/linux/mkl'
mkl_prefix = ['%s/include'%_mklroot, '%s/lib/intel64'%_mklroot]
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']
silo = True
silo_prefix = prefix
silo_libs = ['siloh5', 'hdf5']
trilinos = True
trilinos_prefix = prefix
cppunit_prefix = prefix

env_export = ['INTEL_LICENSE_FILE']
tools_names = [('intelc',{'topdir':'/sw/Licensed/intel/compilers_and_libraries_2016.2.181'})]

