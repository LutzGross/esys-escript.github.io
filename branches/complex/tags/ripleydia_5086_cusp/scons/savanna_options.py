
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

try:
    import os
    if not 'escript/dev-deps' in os.environ['LOADEDMODULES'].split(':'):
        print("WARNING: The escript/dev-deps module does not appear to be loaded!")
except:
    pass

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 201

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
#prefix = '/usr/local'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
#build_dir = 'build'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
#cxx = 'g++'

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
cxx_extra = '-DCORE_ID1 -sox -I/sw/libs/numpy/x86_64/icc-14/1.8-py27_omp/lib/python2.7/site-packages/numpy/core/include'

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = '-shared-intel -L/sw/libs/hdf5/1.8.12-serial/lib -ipo-jobs4'
ld_extra += ' -wd11021 '  #silence icpc warnings about symbols ipo can't see

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
mpi = 'INTELMPI'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
mpi_prefix = '/sw/sdev/intel/impi/4.1.3.049/intel64'

# MPI libraries to link against
#mpi_libs = ['mpi']

# Prefix or paths to boost-python headers and libraries. See note above.
boost_prefix = '/sw/libs/boost/1.55.0-py2.7'

# boost-python library/libraries to link against
boost_libs = ['boost_python']

# Prefix or paths to CppUnit headers and libraries. See note above.
cppunit_prefix = '/sw/apps/cppunit/x86_64/gcc-4.3.2/cppunit-1.12.1'

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
netcdf_prefix = '/sw/libs/netcdf/4.1.3'

# netCDF library/libraries to link against
netcdf_libs = ['netcdf_c++', 'netcdf', 'hdf5']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
parmetis_prefix = '/sw/libs/parmetis/x86_64/icc-13/parmetis-4.0.2'

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
mkl_prefix = ['/sw/sdev/intel/composer_xe_2013_sp1.2.144/mkl/include', '/sw/sdev/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64']

# MKL library/libraries to link against
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']

# Whether to use UMFPACK (requires AMD and BLAS)
# DEFAULT: False
umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
umfpack_prefix = '/sw/libs/umfpack/x86_64/icc-13/umfpack-5.6.1'

# UMFPACK library/libraries to link against
umfpack_libs = ['umfpack', 'amd', 'suitesparseconfig']

# Whether to use BoomerAMG (requires MPI)
# DEFAULT: False
#boomeramg = True

# Prefix or paths to BoomerAMG headers and libraries. See note above.
boomeramg_prefix = '/sw/libs/hypre/x86_64/gcc-4.3.2/hypre-2.0.0'

# BoomerAMG library/libraries to link against
boomeramg_libs = ['HYPRE']
#boomeramg_libs = ['HYPRE_IJ_mv', 'HYPRE_krylov', 'HYPRE_parcsr_ls']

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
# DEFAULT: 'none' (do not use LAPACK)
lapack = 'mkl'

# Prefix or paths to LAPACK headers and libraries. See note above.
lapack_prefix = mkl_prefix

# LAPACK library/libraries to link against
lapack_libs = ['mkl_core']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
silo = True

# Prefix or paths to SILO headers and libraries. See note above.
silo_prefix = '/sw/libs/silo/4.9.1'

# SILO library/libraries to link against
silo_libs = ['siloh5', 'hdf5']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
visit = False

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
visit_prefix = '/sw/apps/visit/2.7.0/linux-x86_64/libsim/V2'

# Sim2 library/libraries to link against
#visit_libs = ['simV2']

# Build dynamic libraries only
#DEFAULT: False
build_shared = True


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

# Use intel's VSL library for random data
# DEFAULT: False
#vsl_random = True

# Extra libraries to link with
#sys_libs = []

# Additional environmental variables to export to the tools
env_export = ['INTEL_LICENSE_FILE']

tools_names = [('intelc',{'topdir':'/sw/sdev/intel/composer_xe_2013_sp1.2.144'})]


#iknowwhatimdoing = False

#forcelazy = 'leave_alone'

#forcecollres = 'leave_alone'

