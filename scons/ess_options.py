# locations of libs etc used by mkl
mkl_path = '/opt/intel/mkl80.019/include'
mkl_lib_path = '/opt/intel/mkl80.019/lib/64'
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
python_path = '/raid2/tools/python-2.3.4/include/python2.3'
python_lib_path = '/raid2/tools/python-2.3.4/lib'
python_lib = 'python2.3'


# locations of libraries for boost
boost_path = '/raid2/tools/boost/include/boost-1_31'
boost_lib_path = '/raid2/tools/boost/lib'
boost_lib = 'boost_python-intel-d-1_31'

# names of libraries for python and boost

# names of c and c++ compilers to use
cc = 'icc'
cxx = 'icc'

# c flags to use
cc_flags  = '-O3 -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -c99 -ansi_alias -no-gcc -fpic -w1'
cc_flags_debug  = '-g -O0 -openmp -openmp_report0 -tpp2 -c99 -ansi_alias -no-gcc -fpic -w1'

# c++ flags to use
cxx_flags = '-O3 -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -ansi -ansi_alias -no-gcc -fpic -w1'
cxx_flags_debug = '-g -O0 -openmp -openmp_report0 -tpp2 -ansi -ansi_alias -no-gcc -fpic -w1 -DDOASSERT -DDOPROF'

# static library archiver flags to use
ar_flags = 'crus'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
