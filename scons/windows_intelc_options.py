
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
cc_flags  = '/EHsc /GR /MD /Qvc7.1'
# Same this does not work. cc_optim  = '/fast /Oi /W3 /Qssp /Qinline-factor-'
cc_optim  = '/Ox /QxP /Qprec-div- /Qssp /Qinline-factor- /Qinline-min-size=0 /Qunroll '
cc_debug  = '/Od /RTCcsu /Zi /Y- /debug:all /Qtrapuv'
omp_optim  = '/Qvec-report0 /Qopenmp /Qopenmp-report0 /Qparallel'
omp_debug  = '/Qvec-report3 /Qopenmp /Qopenmp-report2 /Qparallel'
omp_libs = ['C:\Program Files\Intel\Compiler\C++\9.1\IA32\Lib\libguide40']

tools_names = ['intelc']
