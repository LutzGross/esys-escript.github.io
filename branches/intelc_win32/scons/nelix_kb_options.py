# locations of include files for python
python_path = 'C:/python23/include'
python_lib_path = 'C:/python23/libs'
python_lib = 'python23'

# locations of libraries for boost
boost_path = 'E:/woo409/development/boost_1_33'
boost_lib_path = 'E:/woo409/development/boost_1_33/windows_binary/lib'
boost_lib = 'boost_python-vc71-mt-gd'

cc_defines = ['_USE_MATH_DEFINES', ]
# c flags to use
# 1563 - taking adress of a temporary
# 811 - exception specification for implicitly declared virtual function (destructor usually) incompatible with that of override
# 161 - openmp pargmas are unknown when not compiling with openmp
cc_flags  = '/GR /EHsc /MD /Qc99 /Qwd161 /G7 /O3'
cc_flags_debug  = '/Od /MDd /RTC1 /GR /EHsc /Qc99 /Qwd161'
#cc_flags  = '/GR /EHsc /MD /Qc99 /Qopenmp /Qopenmp-report0 /G7'
#cc_flags_debug  = '/Od /MDd /RTC1 /GR /EHsc /Qc99 /Qopenmp /Qopenmp-report0'

# c++ flags to use
cxx_flags = ''
cxx_flags_debug = ''
# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
sys_libs = []
