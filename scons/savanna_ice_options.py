
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

from savanna_options import *

build_dir = 'buildice'

cxx_extra = '-ipo -sox -I/sw/pymodules/2.7/scipy-0.15.1-ice/lib/python2.7/site-packages/numpy/core/include'

cuda = False
werror = True
