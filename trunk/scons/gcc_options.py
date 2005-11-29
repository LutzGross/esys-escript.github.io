# locations of include files for python and boost
python_path = '/usr/include/python2.3'
boost_path = '/usr/include'

# locations of libraries for python and boost
python_lib_path = '/usr/lib'
boost_lib_path = '/usr/lib'

# names of libraries for python and boost
python_lib = 'python2.3'
boost_lib = 'boost_python'

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

# solver libraries to link with
solver_libs = []
