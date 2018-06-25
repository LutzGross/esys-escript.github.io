
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

# This is a template configuration file for escript on OS X homebrew 11.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
osx_dependency_fix = True
cc_flags = "-pedantic -Wall -fPIC -Wno-unknown-pragmas -Wno-sign-compare -Wno-system-headers -Wno-long-long -Wno-strict-aliasing"
cxx_extra = 'std=c++11 -Wno-c99-extensions'
#mpi = 'OPENMPI'
mpi_prefix = '/usr/local'
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
boost_prefix = '/usr/local'
cppunit_prefix = '/usr/local'
netcdf = True
netcdf_prefix = '/usr/local'
netcdf_libs = ['netcdf_c++', 'netcdf']
werror = False
