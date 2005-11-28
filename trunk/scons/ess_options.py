python_path = '/raid2/tools/python-2.3.4/include/python2.3'
boost_path = '/raid2/tools/boost/include/boost-1_31'

python_lib_path = '/raid2/tools/python-2.3.4/lib'
boost_lib_path = '/raid2/tools/boost/lib'

python_lib = python2.3
boost_lib = boost_python-intel-d-1_31

cc = 'icc'
cxx = 'icc'

cc_flags  = '-O3 -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -c99 -ansi_alias -no-gcc -fpic -w1'
cc_flags_debug  = '-g -O0 -openmp -openmp_report0 -tpp2 -c99 -ansi_alias -no-gcc -fpic -w1'

cxx_flags = '-O3 -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -ansi -ansi_alias -no-gcc -fpic -w1'
cxx_flags_debug = '-g -O0 -openmp -openmp_report0 -tpp2 -ansi -ansi_alias -no-gcc -fpic -w1 -DDOASSERT -DDOPROF'

ar_flags = 'crus'
