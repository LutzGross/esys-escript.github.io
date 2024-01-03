
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

# This is a template configuration file for escript on OS X homebrew.
# Refer to README_FIRST for usage instructions.

import os
escript_opts_version = 203

HOMEBREW_PREFIX = '/opt/homebrew'
pythoncmd = os.path.join(HOMEBREW_PREFIX, 'bin', 'python3')
import subprocess
p=subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1].strip()
revversion=p.stdout.split(' ')[1].split('.')[2].strip()
print(f'Compiled for Python 3.{subversion}.{revversion}.')
pythonlibpath = os.path.join(HOMEBREW_PREFIX,  "Cellar", f"python@3.{subversion}",
                             f"3.{subversion}.{revversion}", "Frameworks", "Python.framework",
                             "Versions", f"3.{subversion}", "lib" )
                             #f"python3.{subversion}", f"config-3.{subversion}-darwin")
pythonincpath = os.path.join(HOMEBREW_PREFIX,  "Cellar", f"python@3.{subversion}",
                             f"3.{subversion}.{revversion}", "Frameworks", "Python.framework",
                             "Versions", f"3.{subversion}", "include", f"python3.{subversion}" )

mpi = 'no'
#mpi_prefix = '/usr/local'
#mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
import glob
GCC=[ os.path.join(HOMEBREW_PREFIX, 'bin', 'g++-13') ]

assert  len(GCC) > 0, "unable to find gcc compiler in "+ os.path.join(HOMEBREW_PREFIX, 'bin')
cxx = GCC[0]
cc = os.path.join(HOMEBREW_PREFIX, 'bin', 'gcc-13')
print(f"c++ compiler is {cxx}")


cppunit_prefix = HOMEBREW_PREFIX

openmp = True

boost_prefix = HOMEBREW_PREFIX
boost_libs = [f'boost_python3{subversion}-mt']
compression_libs = ['boost_iostreams-mt']
#============================================================
netcdf = 4
netcdf_prefix = HOMEBREW_PREFIX
netcdf_libs=['netcdf-cxx4', 'netcdf']
#===========================================================
silo = True
silo_prefix = HOMEBREW_PREFIX
silo_libs = ['siloh5']


lapack =False
lapack_prefix = os.path.join(HOMEBREW_PREFIX, 'lapack' )

umfpack = True
umfpack_prefix = HOMEBREW_PREFIX

build_trilinos = True

#cxx_extra = ''
#ld_extra = '-v'
#cxx = "/opt/homebrew/bin/g++-12"
#cxx = "/usr/bin/clang++"

#cxx = "/opt/local/bin/clang++-mp-14"
#cxx_extra = '-std=c++11'

# LDFLAGS = "-L/opt/homebrew/opt/lapack/lib"
# CPPFLAGS = "-I/opt/homebrew/opt/lapack/include"


#-------


