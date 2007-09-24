
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# flag the MPI settings
# useMPI = 'yes' 
# trilinos_path="/home/Work/trilinos-6/include"
# trilinos_lib_path="/home/Work/trilinos-6/lib"
# trilinos_libs=["aztecoo", "teuchos", "epetra"]


# TODO: Variables named *_path should be *_include

# locations of libs etc used by mkl
# mkl_path = '/opt/intel/mkl80.019/include'
# mkl_lib_path ='/opt/intel/mkl80.019/lib/64'
# mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']

mpi_path = '/home/Work/mpich2-1.0.5p4/include'
mpi_lib_path = '/home/Work/mpich2-1.0.5p4/lib'
mpi_libs = ['mpich', 'rt']
mpi_run = 'mpirun -np 1'

# locations of libs etc used by SCSL
# scsl_path = '/usr/include'
# scsl_lib_path = '/usr/lib'
# scsl_libs = ['scs_mp']
# scsl_libs_MPI = [ 'scs', 'mpi' ]


# locations of include files for python
# python_path = '/data/raid2/toolspp4/python/2.4.1/gcc-3.3.6/include/python2.4'
# python_lib_path = '/data/raid2/toolspp4/python/2.4.1/gcc-3.3.6/lib'
# python_lib = 'python2.4'

# locations of libraries for boost
# boost_path = '/data/raid2/toolspp4/boost/1.33.1/python-2.4.1/gcc-3.3.6/include'
# boost_lib_path = '/data/raid2/toolspp4/boost/1.33.1/python-2.4.1/gcc-3.3.6/lib'
# boost_lib = 'boost_python-mt'

# locations of doc building executables
doxygen_path = '/usr/bin'
epydoc_path = '/usr/bin'

# locations of netcdf
useNetCDF="yes"
netCDF_path = "/home/Work/netcdf-3.6.1/include"
netCDF_lib_path = "/home/Work/netcdf-3.6.1/lib"
netCDF_libs = [ 'netcdf_c++', 'netcdf']

# locations of PAPI
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]

# Comment all this stuff
# c flags to use
# cc_flags  = '-O -fPIC'
# cc_flags_debug  = '-g -fPIC -DTRILINOS -I/home/Work/trilinos-6/include'

# c++ flags to use
# cxx_flags = '-fPIC'
# cxx_flags_debug = '-DDOASSERT -UDOPROF -fPIC -DTRILINOS -I/home/Work/trilinos-6/include'
  	
# c and c++ flags for MPI compilation
# c flags to use
# cc_flags_MPI  = '-O -DPASO_MPI -fPIC'
# cc_flags_debug_MPI  = '-g -DPASO_MPI -fPIC -DTRILINOS -I/home/Work/trilinos-6/include'

# c++ flags to use
# cxx_flags_MPI = '-fPIC'
# cxx_flags_debug_MPI = '-DDOASSERT -UDOPROF -fPIC -DTRILINOS -I/home/Work/trilinos-6/include'

# system specific libraries to link with
# sys_libs = []
