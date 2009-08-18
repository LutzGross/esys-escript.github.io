
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


# Configuration for rodinia (64-bit desktop system running Ubuntu Linux)

import os
ESCRIPT_ROOT		= '/home/Work/InstallArea'

# If you cannot use the default compiler flags set in SConstruct, then change them here
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
# cc_flags		= '-ansi'
# cc_optim		= '-O2'
# cc_debug		= '-g'
# omp_optim		= '-openmp'
# omp_debug		= '-openmp'
# omp_libs		= []

# usedebug		= 'yes'

# Use the default C/C++ flags but add something only for this host:
# cc_extra		= ''
# ld_extra		= ''

# Be picky about errors
usepedantic		= 'yes'

# Extra libraries
# sys_libs		= []

# Python libraries
python_path		= ESCRIPT_ROOT + '/python-2.4.4/include/python2.4'
python_lib_path		= ESCRIPT_ROOT + '/python-2.4.4/lib'
python_libs		= 'python2.4'
python_cmd		= 'python'

# Boost libraries
boost_path		= ESCRIPT_ROOT + '/boost_1_35_0/include/boost-1_35'
boost_lib_path		= ESCRIPT_ROOT + '/boost_1_35_0/lib'
boost_libs		= ['boost_python']

# Specify whether or not to use VTK
usevtk			= 'yes'

# NetCDF
usenetcdf		= 'yes'
netCDF_path		= ESCRIPT_ROOT + '/netcdf-3.6.2/include'
netCDF_lib_path		= ESCRIPT_ROOT + '/netcdf-3.6.2/lib'
netCDF_libs		= ['netcdf_c++', 'netcdf']

# OpenMP
useopenmp		= 'no'

# MPICH2 (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
usempi			= 'no'
mpi_path		= ESCRIPT_ROOT + '/mpich2-1.0.7/include'
mpi_lib_path		= ESCRIPT_ROOT + '/mpich2-1.0.7/lib'
mpi_libs		= ['mpich', 'rt']
mpi_flavour		= 'MPICH'

# MPICH2 for jumpshot (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# usempi		= 'no'
# mpi_path		= ESCRIPT_ROOT + '/mpich2-1.0.7/include'
# mpi_lib_path		= ESCRIPT_ROOT + '/mpich2-1.0.7/lib'
# mpi_libs		= ['lmpe', 'mpe', 'mpich', 'rt']
# mpi_flavour		= 'MPICH'

# ParMETIS (for use with MPI)
useparmetis		= 'yes'
parmetis_path		= ESCRIPT_ROOT + '/parmetis-3.1/include'
parmetis_lib_path	= ESCRIPT_ROOT + '/parmetis-3.1/lib'
parmetis_libs		= ['parmetis', 'metis']

# PAPI
# usepapi		= 'no'
# papi_path		= '/usr/include'
# papi_lib_path		= '/usr/lib'
# papi_libs		= ['papi']
# papi_instrument_solver= 'no'

