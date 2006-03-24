
#          Copyright 2006 by ACcESS MNRF                   
#                                                          
#              http://www.access.edu.au                    
#       Primary Business: Queensland, Australia            
#  Licensed under the Open Software License version 3.0    
#     http://www.opensource.org/licenses/osl-3.0.php       
#                                                          


# locations of libs etc used by mkl
mkl_path = ''
mkl_lib_path = ''
mkl_libs = []

# locations of libs etc used by SCSL
scsl_path = ''
scsl_lib_path = ''
scsl_libs = []

# locations of libs etc used by UMFPACK
umfpack_path = ''
umfpack_lib_path = ''
umfpack_libs = []

# locations of include files for python
python_path = '/opt/tools/python-2.3.5-gcc3.4/include/python2.3'
python_lib_path = '/opt/tools/python-2.3.5-gcc3.4/lib'
python_lib = 'python2.3'

# locations of libraries for boost
boost_path = '/opt/tools/boost_1_31-gcc3.4/include/boost-1_31'
boost_lib_path = '/opt/tools/boost_1_31-gcc3.4/lib'
boost_lib = 'boost_python-gcc-1_31'

# locations of doc building executables
doxygen_path = ''
epydoc_path = ''

# names of c and c++ compilers to use
cc = 'gcc'
cxx = 'g++'

# c flags to use
cc_flags  = '-O3 -std=c99 -fpic --no-warn -W -Wno-unknown-pragmas'
cc_flags_debug  = '-g -O0 -std=c99 -fpic --no-warn -W -Wno-unknown-pragmas'

# c++ flags to use
cxx_flags  = '-O3 -ansi -fpic --no-warn -W -Wno-unknown-pragmas'
cxx_flags_debug  = '-g -O0 -ansi -fpic --no-warn -W -Wno-unknown-pragmas -DDOASSERT -DDOPROF'

# static library archiver flags to use
ar_flags = 'crus'

# system specific libraries to link with
sys_libs = []
