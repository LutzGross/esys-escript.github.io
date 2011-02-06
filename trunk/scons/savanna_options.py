
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

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 200

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
#prefix = '/usr/local'

# Top-level directory for intermediate build and test files.
# DEFAULT: 'build'
#build_dir = 'build'

# C compiler command name or full path.
# DEFAULT: auto-detected
#cc = 'gcc'

# C++ compiler command name or full path.
# DEFAULT: auto-detected
#cxx = 'g++'

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
cc_extra = '-sox' # embed compiler info in binaries

# Additional flags to add to the C++ compiler only
# DEFAULT: '' (empty)
cxx_extra = '-sox' # embed compiler info in binaries

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = '-shared-intel' # fix warning about feupdate in icc v10

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
omp_flags = '-openmp -openmp-report2'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
#omp_ldflags = '-fopenmp'

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
mpi = 'INTELMPI'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
mpi_prefix = ['/sw/sdev/intel/impi/4.0.0.027/include64', '/sw/sdev/intel/impi/4.0.0.027/lib64']

# MPI libraries to link against
mpi_libs = ['mpi']

# Prefix or paths to boost-python headers and libraries. See note above.
boost_prefix = ['/sw/libs/boost/x86_64/gcc-4.1.2/python-2.6.2/boost_1_39_0/include/boost-1_39', '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.6.2/boost_1_39_0/lib']

# boost-python library/libraries to link against
boost_libs = ['boost_python-gcc41-mt']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = True

# Prefix or paths to netCDF headers and libraries. See note above.
netcdf_prefix = '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2'

# netCDF library/libraries to link against
#netcdf_libs = ['netcdf_c++', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
parmetis_prefix = '/sw/libs/parmetis/x86_64/icc-10.1.015/intelmpi/parmetis-3.1.1'

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
mkl_prefix = ['/sw/sdev/intel/mkl/10.2.5.035/include', '/sw/sdev/intel/mkl/10.2.5.035/lib/em64t']

# MKL library/libraries to link against
mkl_libs = ['mkl_core', 'mkl_intel_lp64', 'mkl_intel_thread', 'mkl_lapack', 'guide', 'pthread', 'mkl_mc', 'mkl_def']

# Whether to use UMFPACK (requires AMD and BLAS)
# DEFAULT: False
umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
umfpack_prefix = '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2'

# UMFPACK library/libraries to link against
umfpack_libs = ['umfpack', 'amd', 'blas']

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
# DEFAULT: 'none' (do not use LAPACK)
lapack = 'mkl'

# Prefix or paths to LAPACK headers and libraries. See note above.
lapack_prefix = mkl_prefix

# LAPACK library/libraries to link against
lapack_libs = ['mkl_lapack']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
silo = True

# Prefix or paths to SILO headers and libraries. See note above.
silo_prefix = '/sw/libs/silo/x86_64/gcc-4.3.2/silo-4.7.2'

# SILO library/libraries to link against
silo_libs = ['siloh5', 'hdf5']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
#visit = True

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
visit_prefix = '/sw/apps/visit/x86_64/gcc-4.3.2/visit-2.0.2/2.0.2/linux-x86_64/libsim/V2'

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
env_export = ['INTEL_LICENSE_FILE']

# Build a shared esysUtils library
#share_esysutils = True

# Build a shared paso library
#share_paso = True

tools_names = [('intelc',{'topdir':'/sw/sdev/intel/cc/x86_64/10.1.025'})]

#iknowwhatimdoing = False

#forcelazy = 'leave_alone'

#forcecollres = 'leave_alone'

