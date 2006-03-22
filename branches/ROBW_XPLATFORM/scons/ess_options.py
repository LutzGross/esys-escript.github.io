# locations of libs etc used by mkl
mkl_path = '/opt/intel/mkl80.019/include'
mkl_lib_path ='/opt/intel/mkl80.019/lib/64'
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

# locations of doc building executables
doxygen_path = '/raid2/tools/doxygen/1.4.2/gcc-3.3.5/bin'
epydoc_path = '/raid2/tools/epydoc/2.1/python-2.3.4/bin'
epydoc_pythonpath = '/raid2/tools/epydoc/2.1/python-2.3.4/lib/python2.3/site-packages'

# locations of PAPI
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]
papi_path = ''
papi_lib_path = ''
papi_libs = [ ]

# names of c and c++ compilers to use
cc = 'icc'
cxx = 'icc'

# c flags to use
cc_flags  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -openmp -openmp_report0 -fno-alias -c99 -no-gcc -w1 -fpic"
cc_flags_debug  = '-g -O0 -openmp -openmp_report0 -c99 -ansi_alias -no-gcc -w1 -fpic'

# c++ flags to use
cxx_flags = '-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -openmp -openmp_report0 -ansi -ansi_alias -no-gcc -w1 -fpic'
cxx_flags_debug = '-g -O0 -openmp -openmp_report0 -ansi -ansi_alias -no-gcc -w1  -fpic -DDOASSERT -DDOPROF'

# static library archiver flags to use
ar_flags = 'crus'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
