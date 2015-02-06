
##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
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

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 201

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
prefix = '/opt/escript'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
#build_dir = 'build'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
#cxx = 'g++-4.6'

# Flags to use with the C++ compiler. Do not set unless you know
# what you are doing - use cxx_extra to specify additional flags!
# DEFAULT: compiler-dependent
#cc_flags = ''

# Additional compiler (optimization) flags for non-debug builds
# DEFAULT: compiler-dependent
#cc_optim = '-O3'

# Additional compiler flags for debug builds
# DEFAULT: compiler-dependent
#cc_debug = '-g'

# Additional flags to add to the C++ compiler only
# DEFAULT: '' (empty)
cxx_extra = '-I/usr/lib64/python2.6/site-packages/numpy/core/include -DNO_IMPORT_ARRAY'

# Additional flags to add to the linker
# DEFAULT: '' (empty)
#ld_extra = ''

# Whether to treat compiler warnings as errors
# DEFAULT: True
#werror = False

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
#omp_flags = '-fopenmp'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_ldflags = '-fopenmp'

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
mpi = 'OPENMPI'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
mpi_prefix = ['/usr/include/openmpi-x86_64', '/usr/lib64/openmpi/lib']

# MPI libraries to link against
mpi_libs = ['mpi_cxx', 'mpi', 'dl']

# Prefix or paths to boost-python headers and libraries. See note above.
#boost_prefix = '/usr/local'

# boost-python library/libraries to link against
boost_libs = ['libboost_python-mt']

# Prefix or paths to CppUnit headers and libraries. See note above.
#cppunit_prefix = '/usr/local'

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
#netcdf_prefix = '/usr/local'

# netCDF library/libraries to link against
#netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
#parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
#parmetis_prefix = '/usr/local'

# parMETIS library/libraries to link against
#parmetis_libs = ['parmetis', 'metis']

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
#mkl = True

# Prefix or paths to MKL headers and libraries. See note above.
#mkl_prefix = '/usr'

# MKL library/libraries to link against
#mkl_libs = ['mkl_solver', 'mkl_em64t', 'mkl_core', 'guide', 'pthread']

# Whether to use UMFPACK (requires AMD and BLAS)
# DEFAULT: False
umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib64']

# UMFPACK library/libraries to link against
umfpack_libs = ['umfpack', 'blas']

# Whether to use BoomerAMG (requires MPI)
# DEFAULT: False
#boomeramg = True

# Prefix or paths to BoomerAMG headers and libraries. See note above.
#boomeramg_prefix = '/usr/local'

# BoomerAMG library/libraries to link against
#boomeramg_libs = ['HYPRE']

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
# DEFAULT: 'none' (do not use LAPACK)
lapack = 'clapack'

# Prefix or paths to LAPACK headers and libraries. See note above.
lapack_prefix = ['/usr/include', '/usr/lib64/atlas']

# LAPACK library/libraries to link against
lapack_libs = ['clapack']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
silo = True

# Prefix or paths to SILO headers and libraries. See note above.
silo_prefix = ['/usr/local/include', '/usr/local/lib64']

# SILO library/libraries to link against
silo_libs = ['silo']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
#visit = True

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
#visit_prefix = '/opt/visit/2.4.0/linux-intel/libsim/V2'

# Sim2 library/libraries to link against
#visit_libs = ['simV2']

# Build dynamic libraries only
#DEFAULT: False
#build_shared = True


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

# Use intel's VSL library for random data
# DEFAULT: False
#vsl_random = True

# Extra libraries to link with
#sys_libs = []

# Additional environmental variables to export to the tools
#env_export = []

#tools_names = ['default']

#iknowwhatimdoing = False

#forcelazy = 'leave_alone'

#forcecollres = 'leave_alone'

