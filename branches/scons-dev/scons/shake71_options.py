
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# Configuration for shake71 (Linux)

# Python
python_path		= '/usr/include/python2.5'
python_lib_path		= '/usr/lib/python2.5'
python_lib		= 'python2.5'

# Boost
boost_path		= '/usr/include'
boost_lib_path		= '/usr/lib'
boost_lib		= 'boost_python'

# Documentation
doxygen_path		= '/usr/bin'
epydoc_path		= '/usr/bin'

# NetCDF
useNetCDF		= 'yes'
netCDF_path		= '/usr/include'
netCDF_lib_path		= '/usr/lib'
netCDF_libs		= [ 'netcdf_c++', 'netcdf']

# MPI (version: MPICH2)
mpi_path		= '/home/Work/InstallArea/mpich2-1.0.7/include'
mpi_lib_path		= '/home/Work/InstallArea/mpich2-1.0.7/lib'
mpi_libs		= ['mpich', 'rt']
mpi_run			= 'mpirun -np 1'

# ParMETIS
parmetis_path		= '/home/Work/InstallArea/parmetis-3.1/include'
parmetis_lib_path	= '/home/Work/InstallArea/parmetis-3.1/lib'
parmetis_libs		= ['parmetis', 'metis']

# PAPI
# papi_instrument_solver	= 0
# papi_path		= '/sw/.../include'
# papi_lib_path		= '/sw/.../lib'
# papi_libs		= ['papi']

