
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# flag the MPI settings
useMPI = 'no'

# TODO: Variables named *_path should be *_include

python_version="2.4"
python_installation="2.4.4/icc10.0.023"
boost_version="1_33_1"
boost_installation="1.33.1/python2.4.4/icc9.1.051"

#prefix = ARGUMENTS.get('prefix', '/opt/esys-escript/unstable/')
#tools_prefix="/opt/esys-escript/unstable/"

#    get the installation prefix
# locations of libs etc used by mkl
mkl_path = '/opt/intel_mkl/9.1.018/include'
#mkl_lib_path ='/opt/intel/mkl80.019/lib/64'
mkl_lib_path ='/opt/intel_mkl/9.1.018/lib/64'
mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']

# locations of libs etc used by SCSL
scsl_path = '/opt/scsl/1.6.1.0/include'
scsl_lib_path = '/opt/scsl/1.6.1.0/lib'
scsl_libs = ['scs_mp']
scsl_libs_MPI = [ 'scs', 'mpi' ]


# locations of include files for python
# python_exec= '/opt/python/'+python_installation+'/bin/python'
python_exec= 'python'
python_path = '/opt/python/'+python_installation+'/include/python'+python_version
python_lib_path = '/opt/python/'+python_installation+'/lib'
python_lib = 'python'+python_version

# locations of libraries for boost
boost_path = '/opt/boost/'+boost_installation+'/include/boost-'+boost_version
boost_lib_path = '/opt/boost/'+boost_installation+'/lib'
boost_lib = 'boost_python-mt'

# locations of doc building executables
doxygen_path = '/opt/doxygen-1.4.5/bin'
epydoc_path = '/opt/epydoc-2.1/bin'

# locations of netcdf
useNetCDF = 'yes'
netCDF_path = "/opt/netcdf/3.6.2/include"
netCDF_lib_path = "/opt/netcdf/3.6.2/lib"
netCDF_libs = [ 'netcdf_c++', 'netcdf']

# locations of PAPI
papi_instrument_solver = 0
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]

mpi_path = '/usr/include'
mpi_lib_path = '/usr/lib'
mpi_libs = [ 'mpi' ]
mpi_run = 'mpirun -np 1'

omp_flags = '-openmp -openmp_report2 '
omp_flags_debug = '-openmp -openmp_report0'

# c flags to use
cc_flags  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -c99 -w1 -wd161 -fpic -ivdep-parallel"
cc_flags_debug  = '-g -O0 -c99 -w1 -wd161 -fpic'

# c++ flags to use
cxx_flags = '-ansi -wd161 -DMPI_NO_CPPBIND'
cxx_flags_debug = '-ansi -wd161 -DDOASSERT -DDOPROF -DMPI_NO_CPPBIND'

# c and c++ flags for MPI compilation
# c flags to use
cc_flags_MPI  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -c99 -w1 -fpic -wd161 -DPASO_MPI -ivdep-parallel"
cc_flags_debug_MPI  = '-g -O0 -c99 -w1 -fpic -wd161 -DPASO_MPI'

# c++ flags to use
cxx_flags_MPI = '-ansi -wd1563 -wd161 -DMPI_NO_CPPBIND'
cxx_flags_debug_MPI = '-ansi -DDOASSERT -DDOPROF -wd1563 -wd161 -DMPI_NO_CPPBIND'

# system specific libraries to link with
sys_libs = ['guide', 'irc']

