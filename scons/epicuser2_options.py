
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

# This is a template configuration file for escript/finley on Linux.
# Copy this file to <hostname>_options.py, where <hostname> is your machine's
# short hostname, then customize to your needs.

# PREFIXES:
# There are two ways to specify where to find dependent headers and libraries
# (via the <dependency>_prefix):
# 1) If your installation follows the general scheme where headers are located
#    in <prefix>/include[32,64], and libraries in <prefix>/lib[32,64] then
#    it is sufficient to specify this prefix, e.g. boost_prefix='/usr'
# 2) Otherwise provide a list with two elements, where the first one is the
#    include path, and the second the library path, e.g.
#    boost_prefix=['/usr/include/boost1_44', '/usr/lib']
# All <dependency>_prefix settings default to '/usr'

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 200

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
prefix = '/opt/escript/3.2.1'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
build_dir = 'build'

# C compiler command name or full path.
# DEFAULT: auto-detected
cc = 'gcc'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
cxx = 'g++'

# Flags to use with both C and C++ compilers. Do not set unless you know
# what you are doing - use cc_extra to specify additional flags!
# DEFAULT: compiler-dependent
#cc_flags = ''

# Additional compiler (optimization) flags for non-debug builds
# DEFAULT: compiler-dependent
#cc_optim = '-O3 -mmmx -msse'

# Additional compiler flags for debug builds
# DEFAULT: compiler-dependent
#cc_debug = '-g'

# Additional flags to add to the C compiler only
# DEFAULT: '' (empty)
#cc_extra = ''

# Additional flags to add to the C++ compiler only
# DEFAULT: '' (empty)
#cxx_extra = ''

# Additional flags to add to the linker
# DEFAULT: '' (empty)
#ld_extra = '-shared-intel'

# Whether to treat compiler warnings as errors
# DEFAULT: True
werror = False

# Whether to build a debug version
# DEFAULT: False
debug = False

# Set to True to print the full compiler/linker command line
# DEFAULT: False
#verbose = True

# Set to True to add flags that enable OpenMP parallelization
# DEFAULT: False
openmp = False

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
mpi_prefix = ['/opt/mpi/gcc/openmpi/1.4.3/include', '/opt/mpi/gcc/openmpi/1.4.3/lib']

# MPI libraries to link against
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']

# Prefix or paths to python headers and libraries. See note above.
# By default, this is determined using the running python executable.
python_prefix = ['/opt/python/2.6.7/include/python2.6', '/opt/python/2.6.7/lib']

# Prefix or paths to boost-python headers and libraries. See note above.
boost_prefix = ['/opt/boost/1.39.0/include/boost-1_39', '/opt/boost/1.39.0/lib']

# boost-python library/libraries to link against
boost_libs = ['boost_python-gcc41-mt-1_39']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
netcdf_prefix = ['/opt/netcdf/4.0.1/include', '/opt/netcdf/4.0.1/lib']

# netCDF library/libraries to link against
netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
parmetis_prefix = ['/opt/parmetis/3.1.1/include', '/opt/parmetis/3.1.1/lib']

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
mkl_libs = ['mkl_intel_lp64', 'mkl_gnu_thread', 'libmkl_lapack95_lp64', 'mkl_core', 'gomp', 'pthread']

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
#lapack = 'clapack'

# Prefix or paths to LAPACK headers and libraries. See note above.
#lapack_prefix = '/usr/local'

# LAPACK library/libraries to link against
#lapack_libs = ['lapack_atlas']

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

# Whether to enable the deprecated PyVisi interface (requires the VTK python
# modules)
# DEFAULT: False
#pyvisi = True


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

