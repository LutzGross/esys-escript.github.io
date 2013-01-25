#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

# flag the MPI settings
# useMPI = 'yes' 
# trilinos_path="/usr/local/trilinos-6/include"
# trilinos_lib_path="/usr/local/trilinos-6/lib"
# trilinos_libs=["aztecoo", "teuchos", "epetra"]


# TODO: Variables named *_path should be *_include

# locations of libs etc used by mkl
# mkl_path = '/opt/intel/mkl80.019/include'
# mkl_lib_path ='/opt/intel/mkl80.019/lib/64'
# mkl_libs = ['mkl_solver', 'mkl_lapack', 'mkl_ipf']

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
# doxygen_path = '/raid2/tools/doxygen/1.4.2/gcc-3.3.5/bin'
# epydoc_path = '/raid2/tools/epydoc/2.1/python-2.3.4/bin'
# epydoc_pythonpath = '/raid2/tools/epydoc/2.1/python-2.3.4/lib/python2.3/site-packages'

# locations of PAPI
# papi_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/include'
# papi_lib_path = '/data/raid2/toolspp4/papi/3.0.8.1/gcc-3.3.6/lib'
# papi_libs = [ 'papi' ]

# c flags to use
# cc_flags  = '-O -fPIC'
# cc_flags_debug  = '-g -fPIC'

# c++ flags to use
# cxx_flags = '-fPIC'
# cxx_flags_debug = '-DDOASSERT -UDOPROF -fPIC'
# system specific libraries to link with
sys_libs = []
