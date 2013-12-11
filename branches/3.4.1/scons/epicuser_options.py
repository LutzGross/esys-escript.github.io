
##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################
#

#
#   On epic.ivec.org: link this file to epicuser1_options.py and
#                     epicuser2_options.py 
#

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 201

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
#prefix = '/opt/escript/3.2.1'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
build_dir = 'build'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
cxx = 'icpc'
#cxx = 'g++'

# Flags to use with the C++ compiler. Do not set unless you know
# what you are doing - use cxx_extra to specify additional flags!
# DEFAULT: compiler-dependent
#cc_flags = ''

# Additional compiler (optimization) flags for non-debug builds
# DEFAULT: compiler-dependent
#cc_optim = '-O3 -mmmx -msse'

# Additional compiler flags for debug builds
# DEFAULT: compiler-dependent
cc_debug = '-g'

# Additional flags to add to the C++ compiler only
# DEFAULT: '' (empty)
cxx_extra = '-sox -I/opt/python-numpy/1.6.1-python-2.6/lib/python2.6/site-packages/numpy/core/include'
#cxx_extra = '-fopenmp'

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = '-shared-intel -L/opt/hdf5/1.8.6/lib'

# Whether to treat compiler warnings as errors
# DEFAULT: True
werror = False

# Whether to build a debug version
# DEFAULT: False
#debug = True

# Set to True to print the full compiler/linker command line
# DEFAULT: False
verbose = True

# Set to True to add flags that enable OpenMP parallelization
# DEFAULT: False
openmp = True

# Additional compiler flags for OpenMP builds
# DEFAULT: compiler-dependent
omp_flags = '-openmp -openmp-report2'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_ldflags = '-openmp'

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
mpi = 'OPENMPI'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
#mpi_prefix = '/opt/mpi/intel/openmpi/1.4.3-11.1'
mpi_prefix = '/opt/mpi/gcc/openmpi/1.4.3'

# MPI libraries to link against
#mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']

# Prefix or paths to python headers and libraries. See note above.
# By default, this is determined using the running python executable.
#python_prefix = '/opt/python/2.6.7'

# Prefix or paths to boost-python headers and libraries. See note above.
boost_prefix = '/opt/boost/1.46.1-python-2.6'

# boost-python library/libraries to link against
boost_libs = ['boost_python']

# Prefix or paths to CppUnit headers and libraries. See note above.
cppunit_prefix = '/opt/cppunit/1.12.1'

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']


# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
netcdf_prefix = '/opt/netcdf/4.0.1'

# netCDF library/libraries to link against
netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
parmetis_prefix = '/opt/parmetis/3.1.1'

# parMETIS library/libraries to link against
parmetis_libs = ['parmetis', 'metis']

# Whether to use the Intel PAPI (Performance API) library
# DEFAULT: False
#papi = True

# Prefix or paths to PAPI headers and libraries. See note above.
#papi_prefix = '/usr/local'

# PAPI library/libraries to link against
#papi_libs = ['papi']

# Whether to use PAPI to instrument solver iterations
# DEFAULT: False
#papi_instrument_solver = True

# Whether to use Intel MKL (Math Kernel Library)
# DEFAULT: False
mkl = True

# Prefix or paths to MKL headers and libraries. See note above.
mkl_prefix = ['/opt/intel-mkl/10.3.5.220/mkl/include', '/opt/intel-mkl/10.3.5.220/mkl/lib/intel64']

# MKL library/libraries to link against
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']
#mkl_libs = ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core', 'pthread']

# Whether to use UMFPACK (requires AMD and BLAS)
# DEFAULT: False
#umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
#umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']

# UMFPACK library/libraries to link against
#umfpack_libs = ['umfpack']

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
# DEFAULT: 'none' (do not use LAPACK)
lapack = 'mkl'

# Prefix or paths to LAPACK headers and libraries. See note above.
lapack_prefix = mkl_prefix

# LAPACK library/libraries to link against
#lapack_libs = ['mkl_rt']
lapack_libs = ['mkl_core']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
silo = True

# Prefix or paths to SILO headers and libraries. See note above.
silo_prefix = '/opt/silo/4.7.2'

# SILO library/libraries to link against
silo_libs = ['siloh5', 'hdf5']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
visit = True

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
visit_prefix = '/opt/visit/2.3.1/2.3.1/linux-x86_64/libsim/V2'

# Sim2 library/libraries to link against
#visit_libs = ['simV2']

# Whether to use BoomerAMG (requires MPI)
# DEFAULT: False
#boomeramg = True

# Prefix or paths to BoomerAMG headers and libraries. See note above.
#boomeramg_prefix = '/opt/hypre/2.0.0'

# BoomerAMG library/libraries to link against
#boomeramg_libs = ['HYPRE']
#boomeramg_libs = ['HYPRE_IJ_mv', 'HYPRE_krylov', 'HYPRE_parcsr_ls']


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

# Extra libraries to link with
#sys_libs = []

# Additional environmental variables to export to the tools
#env_export = []

# Build a shared esysUtils library
#share_esysutils = True

# Build a shared paso library
#share_paso = True

#tools_names = ['default']

#iknowwhatimdoing = False

#forcelazy = 'leave_alone'

#forcecollres = 'leave_alone'

