
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

# Replace threaded mkl library by sequential one but don't touch other options
from savanna_options import *
openmp = False
mkl_libs = ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core', 'pthread']

