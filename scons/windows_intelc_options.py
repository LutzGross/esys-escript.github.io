
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
cc_flags  = '/FD /EHsc /GR /Qvc7.1 '
cc_optim  = '/O3 /Qip /Qi /MD /W3 /MD'
cc_debug  = '/Od /RTC1 /MDd /ZI /Yd /Y-'
omp_optim  = '/Qvec-report1 /Qopenmp /Qopenmp-report1 /Qparallel /MD /W3'
omp_debug  = '/Qvec-report3 /Qopenmp /Qopenmp-report2 /Qparallel /MD /W3'
omp_libs = '???'
