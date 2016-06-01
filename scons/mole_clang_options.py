
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

escript_opts_version = 203

cxx_extra = '-std=c++11'
boost_prefix = '/opt/local'
boost_libs = ['boost_python-mt']
cppunit_prefix = '/opt/local'
cppunit_libs = ['cppunit']
netcdf = True
netcdf_prefix = '/opt/local'
netcdf_libs = ['netcdf_c++', 'netcdf']
tools_names = ['clang']

