
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# Configuration for shake71 (Linux)

# ParMETIS
parmetis_path		= '/home/Work/Escript.x86_64/parmetis-3.1/include'
parmetis_lib_path	= '/home/Work/Escript.x86_64/parmetis-3.1/lib'
parmetis_lib		= ['parmetis', 'metis']

# Python
python_path		= '/home/Work/Escript.x86_64/python-2.4.4/include/python2.4'
python_lib_path		= '/home/Work/Escript.x86_64/python-2.4.4/lib'
python_lib		= 'python2.4'

# Boost
boost_path		= '/home/Work/Escript.x86_64/boost_1_33/include/boost-1_33'
boost_lib_path		= '/home/Work/Escript.x86_64/boost_1_33/lib'
boost_lib		= 'boost_python-gcc-mt'

# Documentation
doxygen_path		= '/usr/bin'
epydoc_path		= '/usr/bin'

# NetCDF
useNetCDF		= "yes"
netCDF_path		= "/home/Work/Escript.x86_64/netcdf-3.6.2/include"
netCDF_lib_path		= "/home/Work/Escript.x86_64/netcdf-3.6.2/lib"
netCDF_libs		= [ 'netcdf_c++', 'netcdf']

# MPI
mpi_path		= '/home/Work/mpich2-1.0.5p4/include'
mpi_lib_path		= '/home/Work/mpich2-1.0.5p4/lib'
mpi_libs		= ['mpich', 'rt']
mpi_run			= 'mpirun -np 1'

# PAPI
# papi_instrument_solver	= 0
# papi_path		= '/sw/.../include'
# papi_lib_path		= '/sw/.../lib'
# papi_libs		= ['papi']

# C flags (also used by C++)
# cc_flags		= '-O -fPIC'
# cc_flags_debug	= '-g -fPIC -DTRILINOS -I/home/Work/trilinos-6/include'

# C++ flags
# cxx_flags		= '-fPIC'
# cxx_flags_debug	= '-DDOASSERT -UDOPROF -fPIC -DTRILINOS -I/home/Work/trilinos-6/include'

# System-specific libraries to link with
# sys_libs		= []

