
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

from templates.wheezy_options import *

cxx_extra = '-sox'
ld_extra = '-ipo-separate -shared-intel -diag-disable=10397,11021'
werror = False

boost_libs = ['boost_python-py27']

mkl = True
mkl_prefix = ['/opt/intel/composer_xe_2015/mkl/include', '/opt/intel/composer_xe_2015/mkl/lib/intel64']
mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']

lapack = 'mkl'
lapack_prefix = mkl_prefix
lapack_libs = ['mkl_core']

silo = True
silo_libs = ['siloh5', 'hdf5_openmpi']

tools_names = [('intelc', {'topdir':'/opt/intel/composer_xe_2015'})]

