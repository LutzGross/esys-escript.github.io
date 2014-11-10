
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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

# This is a template configuration file for escript on Ubuntu Linux.
# Create a file named <sourcedir>/scons/<hostname>_options.py, where
# <sourcedir> is the escript source directory and <hostname> is your machine's
# short hostname, add the line
# from templates.lucid_options import *
# then customize to your needs.

# PREFIXES:
# There are two ways to specify where to find dependent headers and libraries
# (via the <dependency>_prefix):
# 1) If your installation follows the general scheme where headers are located
#    in <prefix>/include[32,64], and libraries in <prefix>/lib[32,64] then
#    it is sufficient to specify this prefix, e.g. boost_prefix='/usr'
# 2) Otherwise provide a list with two elements, where the first one is the
#    include path, and the second the library path, e.g.
#    boost_prefix=['/usr/include/boost1_48', '/usr/lib']
# All <dependency>_prefix settings default to '/usr'

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 202

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '<sourcedir>' (source directory)
#prefix = '/usr/local'

# Top-level directory for intermediate build and test files.
# DEFAULT: '<sourcedir>/build'
#build_dir = '/tmp/escriptbuild'

# Set to True to print the full compiler/linker command line
# DEFAULT: False
#verbose = True

# C++ compiler command name or full path.
# DEFAULT: auto-detected
#cxx = 'g++'

# Flags to use with the C++ compiler. Do not set unless you know
# what you are doing - use cxx_extra to specify additional flags!
# DEFAULT: compiler-dependent
#cc_flags = ''

# Additional compiler (optimization) flags for non-debug builds
# DEFAULT: compiler-dependent
#cc_optim = '-O3 -march=native'

# Additional compiler flags for debug builds
# DEFAULT: compiler-dependent
#cc_debug = '-g3 -fno-omit-frame-pointer -D_GLIBCXX_DEBUG'

# Additional flags to add to the C++ compiler
# DEFAULT: '' (empty)
#cxx_extra = '-Wextra -Wno-unused-parameter'

# Additional flags to add to the linker
# DEFAULT: '' (empty)
#ld_extra = ''

# Path to CUDA compiler [new in 202]
# DEFAULT: auto-detected
#nvcc = '/usr/local/bin/nvcc'

# Flags for CUDA compiler [new in 202]
# DEFAULT: '' (empty)
#nvccflags = '-arch=sm_30 -DBOOST_NOINLINE="__attribute__((noinline))"'

# Whether to treat compiler warnings as errors
# DEFAULT: True
#werror = False

# Whether to build a debug version (applying cc_debug flags)
# DEFAULT: False
#debug = True

# Set to True to add flags that enable OpenMP parallelization
# DEFAULT: False
openmp = True

# Additional compiler flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_flags = '-fopenmp'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_ldflags = '-fopenmp'

# Prefix or paths to boost headers and libraries. See note above.
#boost_prefix = '/usr/local'

# boost-python library/libraries to link against
boost_libs = ['boost_python-mt-py26']

# Prefix or paths to CppUnit headers and libraries. See note above.
# Only required for C++ unit tests.
#cppunit_prefix = '/usr/local'

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
#mpi = 'OPENMPI'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
mpi_prefix = '/usr/lib/openmpi'

# MPI libraries to link against
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']

# Whether to add support for GPU-based ripley system matrix (requires nvcc
# and thrust headers) [new in 202]
# DEFAULT: False
#cuda = True

# Prefix or paths to NVidia thrust installation. See note above. [new in 202]
#thrust_prefix = '/usr/local'

# Whether to use the netCDF library for dump file support and netCDF-based
# downunder data import
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
#netcdf_prefix = ['/usr/include/netcdf-3', '/usr/lib']

# netCDF library/libraries to link against
#netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
#parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
#parmetis_prefix = '/usr/local'

# parMETIS library/libraries to link against
#parmetis_libs = ['parmetis', 'metis']

# Whether to add support for the Intel MKL (Math Kernel Library) direct solver
# DEFAULT: False
#mkl = True

# Prefix or paths to MKL headers and libraries. See note above.
#mkl_prefix = ['/opt/intel/composer_xe_2015/mkl/include', '/opt/intel/composer_xe_2015/mkl/lib/intel64']

# MKL library/libraries to link against
#mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']

# Whether to add support for the UMFPACK direct solver (requires AMD and BLAS)
# DEFAULT: False
#umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']

# UMFPACK library/libraries to link against
umfpack_libs = ['umfpack', 'blas', 'amd']

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
#lapack = 'clapack'

# Prefix or paths to LAPACK headers and libraries. See note above.
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']

# LAPACK library/libraries to link against
lapack_libs = ['lapack_atlas']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
#silo = True

# Prefix or paths to SILO headers and libraries. See note above.
#silo_prefix = '/usr/local'

# SILO library/libraries to link against
#silo_libs = ['siloh5', 'hdf5']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
#visit = True

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
#visit_prefix = '/opt/visit/2.1.0/linux-intel/libsim/V2'

# Sim2 library/libraries to link against
#visit_libs = ['simV2']

# List of domain families to build [new in 202]
# DEFAULT: 'all' (i.e. dudley, finley, ripley, speckley)
#domains = ['finley', 'ripley']


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

# Compiler flags for some optimisations in dudley
dudley_assemble_flags = '-funroll-loops'

# enables code that is non-standard
#iknowwhatimdoing = True

# compiler toolset to use
#tools_names = ['intelc']

# Additional environmental variables to export to the tools
#env_export = []

# For testing use only, sets the default value for autolazy
# DEFAULT: 'leave_alone'
#forcelazy = 'on'

# For testing use only, sets the default value for force resolving collective
# operations
# DEFAULT: 'leave_alone'
#forcecollres = 'on'

# Whether to create dynamic libraries for esysUtils and paso
# DEFAULT: False
build_shared = True

# Extra libraries to link with
#sys_libs = []

# Python executable to use for compiling. Must be compatible with the
# boost python library
# DEFAULT: auto-detected (interpreter executing scons)
#pythoncmd = '/usr/bin/python3'

# Whether this is a Python 3 build
# DEFAULT: False
#usepython3 = True

# Name of the python library
# DEFAULT: auto-detected for python 2.x
#pythonlibname = 'python3.4m'

# Path to Python include files
# DEFAULT: auto-detected for python 2.x
#pythonincpath = '/usr/include/python3.4'

# Whether to map index_t to long (for very large matrices) [new in 202]
# DEFAULT: False
#longindices = True

# Enable reading compressed binary grids in ripley? (requires boost iostreams)
# DEFAULT: True
#compressed_files = False

# Compression libraries to link with
# DEFAULT: 'boost_iostreams'
#compression_libs = 'boost_iostreams-mt'

# Whether to use the PAPI (Performance API) library
# DEFAULT: False
#papi = True

# Prefix or paths to PAPI headers and libraries. See note above.
#papi_prefix = '/usr/local'

# PAPI library/libraries to link against
#papi_libs = ['papi']

# Whether to use PAPI to instrument solver iterations
# DEFAULT: False
#papi_instrument_solver = True

