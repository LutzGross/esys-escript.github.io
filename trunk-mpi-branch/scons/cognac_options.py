#          Copyright 2006 by ACcESS MNRF                   
#                                                          
#              http://www.access.edu.au                    
#       Primary Business: Queensland, Australia            
#  Licensed under the Open Software License version 3.0    
#     http://www.opensource.org/licenses/osl-3.0.php       
#                                                          


# locations of libs etc used by mkl
mkl_path = '/opt/intel_mkl/8.0.19/include'
mkl_lib_path = '/opt/intel_mkl/8.0.19/lib/64'
mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']
mpi_run = 'mpirun -np 1'

# locations of libs etc used by SCSL
scsl_path = '/opt/scsl/1.6.1.0/include'
scsl_lib_path = '/opt/scsl/1.6.1.0/lib'
scsl_libs = ['scs_mp']

# locations of include files for python
python_path = '/usr/include/python2.3'
python_lib_path = '/usr/lib'
python_lib = 'python2.3'

# locations of libraries for boost
boost_path = '/opt/boost/python2.3/1.31.0/include'
boost_lib_path = '/opt/boost/python2.3/1.31.0/lib'
boost_lib = 'boost_python-il-mt-1_31'

# locations of doc building executables
doxygen_path = '/opt/doxygen-1.4.5/bin'
epydoc_path = '/opt/epydoc-2.1/bin'
epydoc_pythonpath = '/usr/lib/python2.3/site-packages'

# locations of netcdf
netCDF_path = "/opt/netcdf-3.6.0-p1/intel-9.0/include"
netCDF_lib_path = "/opt/netcdf-3.6.0-p1/intel-9.0/lib"
netCDF_libs_cxx = [ 'netcdf_c++', 'netcdf']

# c flags to use
#cc_flags  = '-O3 -fpic -ip -Ob2 -IPF-fma -ftz -parallel -openmp -mtune=itanium2 -mcpu=itanium2 -c99 -IPF-fltacc -IPF-fp-speculationsafe -ipo -fno-alias'
cc_flags  = '-O3 -fpic -ip -Ob2 -IPF-fma -ftz -parallel -openmp -mtune=itanium2 -mcpu=itanium2 -c99 -IPF-fltacc -IPF-fp-speculationsafe -fno-alias'
cc_flags_debug  = '-g -O0 -fpic -openmp -parallel -c99 -w1'

# c++ flags to use - only need to list the additional ones compared with cc_flags
cxx_flags = '-ansi'
cxx_flags_debug = '-ansi -DDOASSERT -DDOPROF'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
