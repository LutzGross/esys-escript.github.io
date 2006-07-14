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

# locations of libs etc used by SCSL
scsl_path = '/opt/scsl/1.6.1.0/include'
scsl_lib_path = '/opt/scsl/1.6.1.0/lib'
scsl_libs = ['scs_mp']

# locations of include files for python
python_path = '/usr/include/python2.3'
python_lib = 'python2.3'

# locations of libraries for boost
boost_path = '/opt/boost/python2.3/1.33.1/include'
boost_lib_path = '/opt/boost/python2.3/1.33.1/lib'
boost_lib = 'boost_python-gcc-mt-1_33_1'

# c flags to use
cc_flags  = '-O3 -fpic -ip -Ob2 -IPF-fma -ftz -parallel -openmp -mtune=itanium2 -mcpu=itanium2 -c99 -IPF-fltacc -IPF-fp-speculationsafe -ipo -fno-alias'
cc_flags_debug  = '-g -O0 -fpic -openmp -parallel -c99 -w1'

# c++ flags to use - only need to list the additional ones compared with cc_flags
cxx_flags = '-ansi'
cxx_flags_debug = '-ansi -DDOASSERT -DDOPROF'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
