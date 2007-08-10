# locations of include files for python
python_path = 'C:/python23/include'
python_lib_path = 'C:/python23/libs'
python_lib = 'python23'

# locations of libraries for boost
boost_path = '../boost'
boost_libs_path = '../boost/windows_binary/lib'
boost_libs = 'boost_python-vc71-mt-s-1_31'

# locations of netcdf
useNetCDF = "no"
netCDF_path = "..//netcdf/src/include"
netCDF_lib_path = "../netcdf/src/win32/NET/release"
netCDF_libs_cxx = [ 'netcdf', 'netcdf_cpp' ]

cc_defines = ['_USE_MATH_DEFINES', 'BOOST_NO_INTRINSIC_WCHAR_T', 'DLL_NETCDF' ]
# c flags to use
# 1563 - taking adress of a temporary
# 811 - exception specification for implicitly declared virtual function (destructor usually) incompatible with that of override
# 161 - openmp pargmas are unknown when not compiling with openmp
cc_common_flags = '/FD /EHsc /GR /wd4068 '
cc_flags  = cc_common_flags + '/O2 /Op /MT /W3'

cc_flags_debug  = cc_common_flags + '/Od /RTC1 /MTd /ZI'

# c++ flags to use
cxx_flags = ''
cxx_flags_debug = ''
# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
sys_libs = []
