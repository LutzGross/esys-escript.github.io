
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# Configuration for Savanna

# flag the MPI settings
# useMPI = 'yes'

# TODO: Variables named *_path should be *_include

# locations of libs etc used by mkl
### mkl_path = '/sw/sdev/cmkl/10.0.2.18/include'
### mkl_lib_path = '/sw/sdev/cmkl/10.0.2.18/lib/em64t'
### mkl_libs = ['mkl_solver', 'mkl_lapack']
### mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']		# Library mkl_ipf does not exist anymore

# locations of libs etc used by SCSL
### scsl_path = '/usr/include'
### scsl_lib_path = '/usr/lib'
### scsl_libs = ['scs_mp']
### scsl_libs_MPI = [ 'scs', 'mpi' ]

# Location of ParMETIS library
### parmetis_path = '/data/raid2/toolspp4/parmetis/include'
### parmetis_lib_path = '/data/raid2/toolspp4/parmetis/lib'
### parmetis_lib = ['parmetis', 'metis']

# locations of include files for python
python_path = '/sw/apps/python/x86_64/gcc-4.1.2/python-2.4.4/include/python2.4'
python_lib_path = '/sw/apps/python/x86_64/gcc-4.1.2/python-2.4.4/lib'
python_lib = 'python2.4'

# locations of libraries for boost
boost_path = '/usr/include/boost'
boost_lib_path = '/usr/lib64'
boost_lib = 'boost_python'

# locations of doc building executables
### doxygen_path = '/data/raid2/toolspp4/doxygen/1.4.6/gcc-3.3.6/bin'
### epydoc_path = '/raid2/tools/epydoc/2.1/python-2.3.4/bin'

# locations of netcdf
useNetCDF = 'yes'
netCDF_path = '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/include'
netCDF_lib_path = '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/lib'
netCDF_libs = [ 'netcdf_c++', 'netcdf']

# locations of PAPI
papi_instrument_solver = 0
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]

mpi_path = '/usr/include'
mpi_lib_path = '/usr/lib64'
mpi_libs = [ 'mpi' ]
mpi_run = 'mpirun -np 1'

# Change -fno-alias to -fargument-noalias???
# -ivdep-parallel is i64 only, not Savanna
# Use -ax for vectorization?
# -axN and -axW are identical?
# -axK not for SSE2

# Which manual has information about the compiler options?

omp_flags = '-openmp -openmp_report2 '
omp_flags_debug = '-openmp -openmp_report0'

# c flags to use
cc_flags  = '-O3 -ftz -IPF-ftlacc- -IPF-fma -fno-alias -wd161 -fPIC -DBLOCKTIMER'
cc_flags_debug  = '-g -O0 -wd161 -fPIC -DBLOCKTIMER'

# c++ flags to use
cxx_flags = '-ansi -wd161 -DMPI_NO_CPPBIND -DBLOCKTIMER'
cxx_flags_debug = '-g -ansi -wd161 -DDOASSERT -DDOPROF -DMPI_NO_CPPBIND -fPIC -DBLOCKTIMER'

# c and c++ flags for MPI compilation
# c flags to use
cc_flags_MPI  = '-O3 -ftz -IPF-ftlacc- -IPF-fma -fno-alias -wd161 -fPIC -DPASO_MPI -DBLOCKTIMER'
cc_flags_debug_MPI  = '-g -O0 -wd161 -fPIC -DPASO_MPI -DBLOCKTIMER'

# c++ flags to use
cxx_flags_MPI = '-ansi -wd1563 -wd161 -DMPI_NO_CPPBIND -DBLOCKTIMER'
cxx_flags_debug_MPI = '-ansi -DDOASSERT -DDOPROF -wd1563 -wd161 -DMPI_NO_CPPBIND -DBLOCKTIMER'

# system specific libraries to link with
sys_libs = ['stdc++']

