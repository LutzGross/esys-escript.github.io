
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


from windows_mscv71_options import *

# c flags to use
cc_flags  = '/FD /GR /EHs'
cc_optim  = '/O3 /Oi /Qip /MD /W3 /MD'
cc_debug  = '/Od /RTC1 /MDd /ZI /Y-'
omp_optim  = '/Qvec-report0 /Qopenmp /Qopenmp-report0 /Qparallel /MD /W3'
omp_debug  = '/Qvec-report3 /Qopenmp /Qopenmp-report2 /Qparallel /MD /W3'
omp_libs = ['C:\Program Files\Intel\Compiler\C++\9.1\IA32\Lib\libguide']

win_tools_name = 'intelc'
