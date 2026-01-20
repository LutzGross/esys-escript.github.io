
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

import os

# PREFIXES:
# There are two ways to specify where to find dependent headers and libraries
# (via the <dependency>_prefix):
# 1) If your installation follows the general scheme where headers are located
#    in <prefix>/include[32,64], and libraries in <prefix>/lib[32,64] then
#    it is sufficient to specify this prefix, e.g. boost_prefix='C:/python'
# 2) Otherwise provide a list with two elements, where the first one is the
#    include path, and the second the library path, e.g.
#    boost_prefix=['C:/boost/include/boost1_44', 'C:/boost/lib']
# All <dependency>_prefix settings default to '/usr' so have to be set
# manually on Windows.

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
# This setting is mandatory.
escript_opts_version = 200

# Installation prefix. Files will be installed in subdirectories underneath.
# DEFAULT: '.' (current directory)
#prefix = 'C:/escript'

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
cc_flags = '/EHsc /GR /MD /I"C:/Program Files (x86)/Intel/Composer XE 2011 SP1/compiler/include/" /DCPPUNIT_BUILD_DLL'

# Additional compiler (optimization) flags for non-debug builds
# DEFAULT: compiler-dependent
cc_optim = '/fast /Oi /W3 /Qinline-factor- /Qinline-min-size=0 /Qunroll'
cc_optim = '/fast /Oi /Qunroll'

# Additional compiler flags for debug builds
# DEFAULT: compiler-dependent
#cc_debug = '/Od /RTCcsu /Zi /Y- /debug:all /Qtrapuv'

# Additional flags to add to the C compiler only
# DEFAULT: '' (empty)
#cc_extra = ''

# Additional flags to add to the C++ compiler only
# DEFAULT: '' (empty)
#cxx_extra = ''

# Additional flags to add to the linker
# DEFAULT: '' (empty)
ld_extra = '/libpath:"C:/Program Files (x86)/Intel/Composer XE 2011 SP1/compiler/lib/ia32"'

# Whether to treat compiler warnings as errors
# DEFAULT: True
werror = False

# Whether to build a debug version
# DEFAULT: False
#debug = True

# Set to True to print the full compiler/linker command line
# DEFAULT: False
verbose = False

# Set to True to add flags that enable OpenMP parallelization
# DEFAULT: False
openmp = True

# Additional compiler flags for OpenMP builds
# /Qvec-report0 - Vectorise quietly - http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/cpp/lin/copts/common_options/option_vec_report.htm
# /Qopenmp-report0 - OpenMP quietly (No diagnostics report generated.) - http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/cpp/lin/optaps/common/optaps_perf_openmprep.htm
# DEFAULT: compiler-dependent
omp_flags = '/Qvec-report0 /Qopenmp /Qopenmp-report0 /Qparallel'

# Additional linker flags for OpenMP builds
# DEFAULT: compiler-dependent
omp_ldflags = ''

# Flavour of MPI implementation
# Recognized values: 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
# DEFAULT: 'none' (disable MPI)
mpi = 'MPICH2' #'MPICH2'

# Prefix or paths to MPI headers and libraries. See note above about prefixes.
mpi_prefix = "C:/Program Files (x86)/MPICH2"

# MPI libraries to link against
mpi_libs = ['mpi']

dotdot = os.path.realpath('..')

# Prefix or paths to boost-python headers and libraries. See note above.
boost_prefix = ['C:/projects/boost_1_47_0', 'C:/projects/boost_1_47_0/stage/lib']

# boost-python library/libraries to link against
boost_libs = ['boost_python-iw-mt-1_47']

# Prefix or paths to CppUnit headers and libraries. See note above.
#cppunit_prefix = 'C:/CppUnit'

# CppUnit library/libraries to link against
#cppunit_libs = ['cppunit']

# Whether to use the netCDF library for dump file support
# DEFAULT: False
netcdf = False

# Prefix or paths to netCDF headers and libraries. See note above.
netcdf_prefix = [os.path.join(dotdot, 'netcdf', 'src', 'include'), os.path.join(dotdot, 'netcdf', 'lib')]

# netCDF library/libraries to link against
netcdf_libs = ['netcdf_cpp', 'netcdf']

# Whether to use the parMETIS library (only in conjunction with MPI)
# DEFAULT: False
#parmetis = True

# Prefix or paths to parMETIS headers and libraries. See note above.
#parmetis_prefix = 'C:/parmetis'

# parMETIS library/libraries to link against
#parmetis_libs = ['parmetis', 'metis']

# Whether to use the Intel PAPI (Performance API) library
# DEFAULT: False
#papi = True

# Prefix or paths to PAPI headers and libraries. See note above.
#papi_prefix = 'C:/papi'

# PAPI library/libraries to link against
#papi_libs = ['papi']

# Whether to use PAPI to instrument solver iterations
# DEFAULT: False
#papi_instrument_solver = True

# Whether to use Intel MKL (Math Kernel Library)
# DEFAULT: False
mkl = True

# Prefix or paths to MKL headers and libraries. See note above.
mkl_base_dir = 'C:/Program Files (x86)/Intel/Composer XE 2011 SP1/mkl'
mkl_prefix = [os.path.join(mkl_base_dir, 'include'), os.path.join(mkl_base_dir, 'lib/ia32')]

# MKL library/libraries to link against
mkl_libs = ['mkl_intel_c_dll', 'mkl_intel_thread_dll', 'mkl_core_dll', 'libiomp5md', 'mkl_solver', 'mkl_intel_thread', 'mkl_core']

# Whether to use UMFPACK (requires AMD and BLAS)
# DEFAULT: False
#umfpack = True

# Prefix or paths to UMFPACK headers and libraries. See note above.
#umfpack_prefix = 'C:/umfpack'

# UMFPACK library/libraries to link against
#umfpack_libs = ['umfpack']

# Whether to use BoomerAMG (requires MPI)
# DEFAULT: False
#boomeramg = True

# Prefix or paths to BoomerAMG headers and libraries. See note above.
#boomeramg_prefix = 'C:/boomeramg'

# BoomerAMG library/libraries to link against
#boomeramg_libs = ['HYPRE']

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
# DEFAULT: 'none' (do not use LAPACK)
#lapack = 'clapack'
lapack = 'mkl'

# Prefix or paths to LAPACK headers and libraries. See note above.
#lapack_prefix = 'C:/lapack'
lapack_prefix = mkl_prefix

# LAPACK library/libraries to link against
lapack_libs = ['mkl_lapack95']

# Whether to use LLNL's SILO library for Silo output file support in weipa
# DEFAULT: False
#silo = True

# Prefix or paths to SILO headers and libraries. See note above.
#silo_prefix = 'C:/silo'

# SILO library/libraries to link against
#silo_libs = ['siloh5', 'hdf5']

# Whether to use LLNL's VisIt simulation interface (only version 2 supported)
# DEFAULT: False
#visit = True

# Prefix or paths to VisIt's sim2 headers and libraries. See note above.
#visit_prefix = 'C:/visit/2.1.0/linux-intel/libsim/V2'

# Sim2 library/libraries to link against
#visit_libs = ['simV2']


### ADVANCED OPTIONS ###
# Do not change the following options unless you know what they do

# Use intel's VSL library for random data
# DEFAULT: False
#vsl_random = True

# Extra libraries to link with
# sys_libs = ['C:/Program Files/Microsoft Visual Studio .NET 2003/Vc7/PlatformSDK/Lib/Ws2_32', 'C:\Program Files\Intel\Compiler\C++\9.1\IA32\Lib\libguide40']
sys_libs = ['']

# Additional environmental variables to export to the tools
#env_export = []

# overloading - [('intelc',{'topdir':'/sw/sdev/intel/cc/x86_64/10.1.025'})]
tools_names = [('intelc',{'topdir':'C:/Program Files (x86)/Intel/Composer XE 2011 SP1'})]

#iknowwhatimdoing = False

#forcelazy = 'leave_alone'

#forcecollres = 'leave_alone'

