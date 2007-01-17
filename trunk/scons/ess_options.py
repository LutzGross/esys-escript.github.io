
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# flag the MPI settings
# useMPI = 'yes' 

# TODO: Variables named *_path should be *_include

# locations of libs etc used by mkl
mkl_path = '/opt/intel/mkl80.019/include'
mkl_lib_path ='/opt/intel/mkl80.019/lib/64'
mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']

# locations of libs etc used by SCSL
scsl_path = '/usr/include'
scsl_lib_path = '/usr/lib'
scsl_libs = ['scs_mp']
scsl_libs_MPI = [ 'scs', 'mpi' ]


# locations of include files for python
python_path = '/data/raid2/toolspp4/python/2.4.1/gcc-3.3.6/include/python2.4'
python_lib_path = '/data/raid2/toolspp4/python/2.4.1/gcc-3.3.6/lib'
python_lib = 'python2.4'

# locations of libraries for boost
boost_path = '/data/raid2/toolspp4/boost/1.33.1/python-2.4.1/gcc-3.3.6/include'
boost_lib_path = '/data/raid2/toolspp4/boost/1.33.1/python-2.4.1/gcc-3.3.6/lib'
boost_lib = 'boost_python-mt'

# locations of doc building executables
doxygen_path = '/raid2/tools/doxygen/1.4.2/gcc-3.3.5/bin'
epydoc_path = '/raid2/tools/epydoc/2.1/python-2.3.4/bin'
epydoc_pythonpath = '/raid2/tools/epydoc/2.1/python-2.3.4/lib/python2.3/site-packages'

# locations of PAPI
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]

# c flags to use
cc_flags  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -openmp -openmp_report2 -c99 -w1 -fpic -ivdep-parallel"
cc_flags_debug  = '-g -O0 -openmp -openmp_report0 -c99 -w1 -fpic'

# c++ flags to use
cxx_flags = '-ansi'
cxx_flags_debug = '-ansi -DDOASSERT -DDOPROF'
  	
# c and c++ flags for MPI compilation
# c flags to use
cc_flags_MPI  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -c99 -w1 -fpic -wd161 -DPASO_MPI -ivdep-parallel"
cc_flags_debug_MPI  = '-g -O0 -c99 -w1 -fpic -wd161 -DPASO_MPI'

# c++ flags to use
cxx_flags_MPI = '-ansi -wd1563 -wd161'
cxx_flags_debug_MPI = '-ansi -DDOASSERT -DDOPROF -wd1563 -wd161'

# system specific libraries to link with
sys_libs = ['guide', 'irc']
