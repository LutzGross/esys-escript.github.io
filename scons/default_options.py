python_path = '/usr/include'
boost_path = '/usr/include'

python_lib_path = '/usr/lib'
boost_lib_path = '/usr/lib'

python_lib = 'python2.3'
boost_lib = 'boost_python'

cc = 'gcc'
cxx = 'g++'

cc_flags  = '-O3 -std=c99 -fpic -W -Wall -Wno-unknown-pragmas'
cc_flags_debug  = '-g -O0 -std=c99 -fpic -W -Wall -Wno-unknown-pragmas'

cc_flags  = '-O3 -ansi -fpic -W -Wall -Wno-unknown-pragmas'
cc_flags_debug  = '-g -O0 -ansi -fpic -W -Wall -Wno-unknown-pragmas -DDOASSERT -DDOPROF'

ar_flags = 'crus'
