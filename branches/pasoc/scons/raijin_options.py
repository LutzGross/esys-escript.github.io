
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

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 202

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
#prefix = '/usr/local'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
build_dir = '/short/r17/cha564/escript_build'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
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
#cc_debug = '-g'

# Additional flags to add to the C++ compiler
# DEFAULT: '' (empty)
cxx_extra = '-sox -I/apps/python/2.7.3/lib/python2.7/site-packages/numpy/core/include -I/apps/metis/5.0.2/include -wd981 -pthread'

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = '-shared-intel -pthread -L/apps/metis/5.0.2/lib'

# Whether to treat compiler warnings as errors
# DEFAULT: True
werror = False

# Whether to build a debug version
# DEFAULT: False
# debug = True

# Set to True to print the full compiler/linker command line
# DEFAULT: False
verbose = True

# Set to True to add flags that enable OpenMP parallelization
# DEFAULT: False
openmp = True

# Additional compiler flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_flags = '-openmp -openmp-report=1'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_ldflags = '-fopenmp'

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
mpi = 'OPENMPI'
#mpi = 'INTELMPI'

#icpc -I/apps/openmpi/1.6.5/include -pthread -L/apps/openmpi/1.6.5/lib -lmpi_cxx -lmpi -ldl -lm -lnuma -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ldl

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
mpi_prefix = '/apps/openmpi/1.6.5'
#mpi_prefix = ['/apps/intel-mpi/4.1.1.036/include64', '/apps/intel-mpi/4.1.1.036/lib64']

# MPI libraries to link against
mpi_libs = ['mpi_cxx', 'mpi', 'dl', 'm', 'numa', 'rt', 'nsl', 'util']


# Prefix or paths to boost-python headers and libraries. See note above.
boost_prefix = ['/apps/boost/1.54.0-python27/include', '/apps/boost/1.54.0-python27/lib/Intel']

# boost-python library/libraries to link against
boost_libs = ['boost_python']

# Prefix or paths to CppUnit headers and libraries. See note above.
cppunit_prefix = '/apps/cppunit/1.12.1'

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
netcdf_prefix = '/apps/netcdf/4.1.3'

# netCDF library/libraries to link against
#netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
parmetis_prefix = '/apps/parmetis/4.0.2'

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
mkl = True

# Prefix or paths to MKL headers and libraries. See note above.
mkl_prefix = ['/apps/intel-ct/13.4.183/mkl/include', '/apps/intel-ct/13.4.183/mkl/lib/intel64']

# MKL library/libraries to link against
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']

# Whether to use UMFPACK (requires AMD and BLAS)
# DEFAULT: False
#umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
#umfpack_prefix = '/sw/libs/umfpack/x86_64/icc-13/umfpack-5.6.1'

# UMFPACK library/libraries to link against
#umfpack_libs = ['umfpack', 'amd', 'suitesparseconfig']

# Whether to use BoomerAMG (requires MPI)
# DEFAULT: False
#boomeramg = True

# Prefix or paths to BoomerAMG headers and libraries. See note above.
#boomeramg_prefix = '/sw/libs/hypre/x86_64/gcc-4.3.2/hypre-2.0.0'

# BoomerAMG library/libraries to link against
#boomeramg_libs = ['HYPRE']
#boomeramg_libs = ['HYPRE_IJ_mv', 'HYPRE_krylov', 'HYPRE_parcsr_ls']

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
# DEFAULT: 'none' (do not use LAPACK)
lapack = 'mkl'
#lapack = 'none'

# Prefix or paths to LAPACK headers and libraries. See note above.
lapack_prefix = mkl_prefix

# LAPACK library/libraries to link against
lapack_libs = ['mkl_core']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
silo = True

# Prefix or paths to SILO headers and libraries. See note above.
silo_prefix = '/apps/silo/4.9.1'

# SILO library/libraries to link against
silo_libs = ['silo']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
#visit = False

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
#visit_prefix = '/sw/apps/visit/x86_64/gcc-4.3.2/visit-2.6.0/2.6.0/linux-x86_64/libsim/V2'

# Sim2 library/libraries to link against
#visit_libs = ['simV2']


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

# Use intel's VSL library for random data
# DEFAULT: False
#vsl_random = True

# Extra libraries to link with
#sys_libs = []

# Additional environmental variables to export to the tools
env_export = ['INTEL_LICENSE_FILE']

tools_names = [('intelc',{'topdir':'/apps/intel-ct/13.4.183'})]

#iknowwhatimdoing = False

#forcelazy = 'leave_alone'

#forcecollres = 'leave_alone'

