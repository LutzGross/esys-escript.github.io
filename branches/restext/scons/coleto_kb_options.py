
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


# locations of include files for python
python_path = 'C:/python23/include'
python_lib_path = 'C:/python23/libs'
python_lib = 'python23'

# locations of libraries for boost
boost_path = 'c:/woo409/development/boost_1_33'
boost_lib_path = 'c:/woo409/development/boost_1_33/windows_binary/lib'
boost_lib = 'boost_python-vc71-mt-gd'

cc_defines = ['_USE_MATH_DEFINES', ]
# c flags to use
# 1563 - taking adress of a temporary
# 811 - exception specification for implicitly declared virtual function (destructor usually) incompatible with that of override
# 161 - openmp pargmas are unknown when not compiling with openmp
cc_flags  = '/GR /EHsc /MD /Qc99 /Qopenmp /Qopenmp-report1 /O3 /G7 /Qprec /Qparallel /Qpar-report1 /QxP /QaxP'
#cc_flags  = '/GR /EHsc /MD /Qc99 /O3 /G7 /Qprec /QxP /QaxP'

cc_flags_debug  = '/Od /MDd /RTC1 /GR /EHsc /Qc99 /Qopenmp /Qopenmp-report1 /Qprec'

# c++ flags to use
cxx_flags = ''
cxx_flags_debug = ''
# static library archiver flags to use
#ar_flags = 'crus'

# system specific libraries to link with
sys_libs = []
