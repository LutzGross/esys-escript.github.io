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
scsl_path = ''
scsl_lib_path = ''
scsl_libs = ['scs_mp']

                                                                                                                                                                                                     
# locations of libs etc used by UMFPACK
umfpack_path = ''
umfpack_lib_path = ''
umfpack_libs = []

# locations of include files for python
python_path = '/usr/include/python2.3'
python_lib_path = '/usr/lib/python2.3/lib'
python_lib = 'python2.3'


# locations of libraries for boost
boost_path = '/home/woo409/dev/boost_1_33_1'
boost_lib_path = '/home/woo409/dev/boost_1_33_1/altix_binary/lib'
boost_lib = 'boost_python-il-mt-1_33_1'

# names of libraries for python and boost

# names of c and c++ compilers to use
cc = 'icc'
cxx = 'icpc'

# c flags to use
cc_flags  = '-O3 -fpic -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -c99 -ansi_alias -w1'
cc_flags_debug  = '-g -O0 -fpic -openmp -openmp_report0 -tpp2 -c99 -ansi_alias -w1'

# c++ flags to use
cxx_flags = '-O3 -fpic -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -ansi -ansi_alias -w1'
cxx_flags_debug = '-g -O0 -fpic -openmp -openmp_report0 -tpp2 -ansi -ansi_alias -w1 -DDOASSERT -DDOPROF'

# static library archiver flags to use
ar_flags = 'crus'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
