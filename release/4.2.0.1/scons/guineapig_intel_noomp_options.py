
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

# Replace threaded mkl library by sequential one but don't touch other options
from guineapig_intel_options import *
openmp = False
mkl_libs = ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core', 'pthread']

# prevent icpc from complaining about unrecognized omp pragmas
cxx_extra += ' -diag-disable=3180'

