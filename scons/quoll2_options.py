escript_opts_version=203
pythoncmd='/usr/bin/python3'
pythonlibname='python3.7m'
pythonlibpath='/usr/lib/x86_64-linux-gnu/'
pythonincpath='/usr/include/python3.7m'
boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
boost_libs='boost_python37'
debug=1
cxx_extra+=" -O0 -p -pg "
cxx_extra+=" -Wno-error=unused-variable -Wno-error=unused-but-set-variable"
cxx_extra+=" -Wno-error=comment "
werror=0
paso=0
trilinos=1
trilinos_prefix='/usr/local/trilinos_nompi/include'