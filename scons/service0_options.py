
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# Configuration for Savanna (SGI ICE 8200)

# ParMETIS
# parmetis_path		= '/sw/libs/parmetis/x86_64/gcc-4.1.2/parmetis-3.1/include'
# parmetis_lib_path	= '/sw/libs/parmetis/x86_64/gcc-4.1.2/parmetis-3.1/lib'
# parmetis_lib		= ['parmetis', 'metis']

# Python
python_path		= '/sw/apps/python/x86_64/gcc-4.1.2/python-2.4.4/include/python2.4'
python_lib_path		= '/sw/apps/python/x86_64/gcc-4.1.2/python-2.4.4/lib'
python_lib		= 'python2.4'

# Boost
boost_path		= '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.4.4/boost_1_33/include/boost-1_33'
boost_lib_path		= '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.4.4/boost_1_33/lib'
boost_lib		= 'boost_python-gcc'

# Documentation
# doxygen_path		= '/sw/apps/.../bin'
# epydoc_path		= '/sw/apps/.../bin'

# NetCDF
useNetCDF		= 'yes'
netCDF_path		= '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/include'
netCDF_lib_path		= '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/lib'
netCDF_libs		= ['netcdf_c++', 'netcdf']

# MPI
mpi_path		= '/usr/include'
mpi_lib_path		= '/usr/lib64'
mpi_libs		= ['mpi']
mpi_run			= 'mpirun -np 1'

# PAPI
papi_instrument_solver	= 0
# papi_path		= '/sw/.../include'
# papi_lib_path		= '/sw/.../lib'
# papi_libs		= ['papi']

# MKL
# mkl_path		= '/sw/sdev/cmkl/10.0.2.18/include'
# mkl_lib_path		= '/sw/sdev/cmkl/10.0.2.18/lib/em64t'
# mkl_libs		= ['mkl_solver', 'mkl_lapack']

# OpenMP (comment out to disable OpenMP)
omp_flags		= '-openmp -openmp_report0'
omp_flags_debug		= '-openmp -openmp_report0'

# C flags (also used by C++)
cc_flags		= '-O3 -ansi -fPIC -vec-report0 -ftz -IPF-fltacc- -IPF-fma -fno-alias -DBLOCKTIMER -UPASO_DYNAMIC_SCHEDULING_MVM'
cc_flags_debug		= '-g  -ansi -fPIC'

# C++ flags
cxx_flags		= ''
cxx_flags_debug		= '-DDOASSERT -DDOPROF'	# -D... here is not recognized by scons as dependencies

# System-specific libraries to link with
sys_libs		= ['guide', 'pthread', 'stdc++']

