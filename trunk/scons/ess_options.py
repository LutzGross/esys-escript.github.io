# locations of include files for python and boost
python_path = '/raid2/tools/python-2.3.4/include/python2.3'
boost_path = '/raid2/tools/boost/include/boost-1_31'

# locations of libraries for python and boost
python_lib_path = '/raid2/tools/python-2.3.4/lib'
boost_lib_path = '/raid2/tools/boost/lib'

# names of libraries for python and boost
python_lib = 'python2.3'
boost_lib = 'boost_python-intel-d-1_31'

# names of c and c++ compilers to use
cc = 'icc'
cxx = 'icc'

# c flags to use
cc_flags  = '-O3 -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -c99 -ansi_alias -no-gcc -fpic -w1 -DSCSL'
cc_flags_debug  = '-g -O0 -openmp -openmp_report0 -tpp2 -c99 -ansi_alias -no-gcc -fpic -w1 -DSCSL'

# c++ flags to use
cxx_flags = '-O3 -IPF_fma -ftz -openmp -openmp_report0 -mp1 -tpp2 -ansi -ansi_alias -no-gcc -fpic -w1 -DSCSL'
cxx_flags_debug = '-g -O0 -openmp -openmp_report0 -tpp2 -ansi -ansi_alias -no-gcc -fpic -w1 -DDOASSERT -DDOPROF -DSCSL'

# static library archiver flags to use
ar_flags = 'crus'

# system specific libraries to link with
sys_libs = ['guide', 'irc']

# solver libraries to link with
solver_libs = ['scs_mp']
