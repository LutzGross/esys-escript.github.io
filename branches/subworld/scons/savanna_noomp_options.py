
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

# Replace threaded mkl library by sequential one but don't touch other options
from savanna_options import *
openmp = False
mkl_libs = ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core', 'pthread']

