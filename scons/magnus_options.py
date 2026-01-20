
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 203

import os

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
#prefix = '/usr/local'
#prefix = '/group/geosciences953/escript_directory/escript'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
#build_dir = 'build'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
cxx = 'CC'

# Flags to use with the C++ compiler. Do not set unless you know
# what you are doing - use cc_extra to specify additional flags!
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
#cxx_extra = '-shared -fPIC -h gnu -h nomessage=47:1199:1794:1836:11709'
cxx_extra = '-fPIC'

# Additional flags to add to the linker
# DEFAULT: '' (empty)
#ld_extra = '-shared-intel -L/opt/cray/hdf5/1.8.11/cray/81/lib -ipo-jobs4'
#ld_extra = '-dynamic -L/ivec/cle50/devel/PrgEnv-gnu/5.0.41/hdf5/1.8.12/lib'
ld_extra = '-dynamic'

# launcher, prelaunch, postlaunch: for MPI builds/batch system runs
# the following substitutions are applied to all three:
# %b = executable, %n = number of nodes, %p = number of processes,
# %N = total number of processes, # %t = number of threads,
# %f = name of hostfile, %h = comma-separated list of hosts,
# %e = comma-separated list of environment variables to export
#prelaunch = "EE=$(echo -x %e|sed -e 's/,/ -x /g')"
prelaunch = ""
launcher = "aprun -B %b"
postlaunch = ""

stdlocationisprefix = True

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
#omp_flags = '-homp'
omp_flags = '-fopenmp'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
omp_ldflags = '-fopenmp'

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
mpi = 'MPICH2'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
# on magnus, this variable should be set by the programming environment
mpi_prefix = os.environ['MPICH_DIR']


# MPI libraries to link against. Compiler wrapper takes care of this
#mpi_libs = ['mpich', 'mpichcxx']
mpi_libs = []


# Prefix or paths to boost-python headers and libraries. See note above.
#boost_prefix = '/ivec/cle50/devel/PrgEnv-intel/boost/1.55.0'
#boost_prefix = '/ivec/cle50/devel/PrgEnv-gnu/5.0.41/boost/1.49.0'
#boost_prefix = ['/home/caltinay/boost_1_55_0','/home/caltinay/boost_1_55_0/stage/lib']
boost_prefix = '/group/pawsey0143/software/cle52up04/apps/PrgEnv-gnu/5.2.82/gcc/4.9.2/haswell/boost/1.57.0'
# boost-python library/libraries to link against
boost_libs = ['boost_python']

# Prefix or paths to CppUnit headers and libraries. See note above.
#cppunit_prefix = ''

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
#netcdf_prefix = '/opt/cray/netcdf-hdf5parallel/4.3.0/CRAY/81'
#netcdf_prefix = '/ivec/cle50/devel/PrgEnv-gnu/5.0.41/netcdf/4.1.3'
#netcdf_prefix = '/group/geosciences953/escript_directory/escript'
netcdf_prefix = '/group/pawsey0143/joel/netcdf'
#netcdf_prefix=os.environ["NETCDF_DIR"]

# netCDF library/libraries to link against
netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
#parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
#parmetis_prefix = '/sw/libs/parmetis/x86_64/icc-13/parmetis-4.0.2'

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
#mkl_prefix = ['/opt/intel/composer_xe_2013.5.192/mkl/include', '/opt/intel/composer_xe_2013.5.192/mkl/lib/intel64']

# MKL library/libraries to link against
#mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']

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
#lapack = 'mkl'

# Prefix or paths to LAPACK headers and libraries. See note above.
#lapack_prefix = mkl_prefix

# LAPACK library/libraries to link against
#lapack_libs = ['mkl_core']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
#silo = True

# Prefix or paths to SILO headers and libraries. See note above.
silo_prefix = '/group/geosciences953/escript_directory/escript'

# SILO library/libraries to link against
silo_libs = ['silo']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
#visit = True

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
#visit_prefix = ''

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
# On Cray the compiler wrapper depends on a lot of environment vars
# so we simply export everything
import os
env_export = os.environ.keys()

#tools_names = [('intelc',{'topdir':'/opt/intel/composer_xe_2013.5.192'})]

#iknowwhatimdoing = False

#forcelazy = 'auto'

#forcecollres = 'auto'

