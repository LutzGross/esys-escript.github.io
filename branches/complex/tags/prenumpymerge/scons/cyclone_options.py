
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


# flag the MPI settings
# useMPI = 'yes'

# TODO: Variables named *_path should be *_include

# locations of libs etc used by mkl (cyclone is missile mkl_solver.h)
# mkl_path = '/HPC/apps/Intel/mkl701/lib/64'
# mkl_lib_path ='/HPC/apps/Intel/mkl701/include'
# mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']

# locations of libs etc used by SCSL
scsl_path = '/usr/include'
scsl_lib_path = '/usr/lib'
scsl_libs = ['scs_mp']
scsl_libs_MPI = [ 'scs', 'mpi' ]


# locations of include files for python
python_path = '/home/uqksteub/Escript/Python-2.4.4/include/python2.4'
python_lib_path = '/home/uqksteub/Escript/Python-2.4.4/lib'
python_lib = 'python2.4'

# locations of libraries for boost
boost_path = '/HPC/home/uqksteub/Escript/boost/include'
boost_lib_path = '/HPC/home/uqksteub/Escript/boost/lib'
boost_lib = 'boost_python-mt'

# locations of doc building executables
doxygen_path = '/data/raid2/toolspp4/doxygen/1.4.6/gcc-3.3.6/bin'
epydoc_path = '/raid2/tools/epydoc/2.1/python-2.3.4/bin'

# locations of netcdf
useNetCDF = 'yes'
netCDF_path = '/HPC/apps/ia64/include'
netCDF_lib_path = '/HPC/apps/ia64/lib'
netCDF_libs = [ 'netcdf_c++', 'netcdf']

# locations of PAPI
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]

# c flags to use
cc_flags  = '-O -fPIC -openmp -openmp_report2'
cc_flags_debug  = '-g -fPIC -openmp -openmp_report2'

# c++ flags to use
cxx_flags = '-fPIC -openmp -openmp_report2'
cxx_flags_debug = '-DDOASSERT -UDOPROF -fPIC -openmp -openmp_report2'

# c and c++ flags for MPI compilation
# c flags to use
### cc_flags_MPI  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -fpic -wd161 -DPASO_MPI -ivdep-parallel"
### cc_flags_debug_MPI  = '-g -O0 -fpic -wd161 -DPASO_MPI'

# c++ flags to use
### cxx_flags_MPI = '-ansi -wd1563 -wd161'
### cxx_flags_debug_MPI = '-ansi -DDOASSERT -DDOPROF -wd1563 -wd161'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
