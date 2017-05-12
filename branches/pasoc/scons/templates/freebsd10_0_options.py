
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

# This is a template configuration file for escript on FreeBSD.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
cxx_extra = '-I/usr/local/lib/python2.7/site-packages/numpy/core/include'
openmp = True
boost_prefix = '/usr/local'
boost_libs = ['boost_python']
cppunit_prefix = '/usr/local'
netcdf = True
netcdf_prefix = '/usr/local'
silo = True
silo_prefix = '/usr/local'
silo_libs = ['silo']

