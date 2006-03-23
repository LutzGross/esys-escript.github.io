import sys

# locations of include files for python
python_path = "/usr/include/python%s.%s"%(sys.version_info[0],sys.version_info[1])
python_lib_path = '/usr/lib'
python_lib = "python%s.%s"%(sys.version_info[0],sys.version_info[1])

# locations of libraries for boost
boost_path = '/usr/include'
boost_lib_path = '/usr/lib'
boost_lib = 'boost_python'

# locations of doc building executables
doxygen_path = '/usr/bin'
epydoc_path = '/usr/bin'

# c flags to use
cc_flags  = '-O3 -std=c99 -fpic --no-warn -W -Wno-unknown-pragmas'
cc_flags_debug  = '-g -O0 -std=c99 -fpic --no-warn -W -Wno-unknown-pragmas'

# c++ flags to use
cxx_flags  = '-O3 -ansi -fpic --no-warn -W -Wno-unknown-pragmas'
cxx_flags_debug  = '-g -O0 -ansi -fpic --no-warn -W -Wno-unknown-pragmas -DDOASSERT -DDOPROF'

# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
#sys_libs = []
