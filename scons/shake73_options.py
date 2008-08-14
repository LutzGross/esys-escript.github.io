
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# Configuration for Savanna (SGI ICE 8200 running SUSE Linux)

# If you cannot use the default compiler flags set in SConstruct, then change them here
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
# cc_flags		= '-ansi'
# cc_optim		= '-O2'
# cc_debug		= '-g'
# omp_optim		= '-openmp'
# omp_debug		= '-openmp'
# omp_libs		= []

# Use the default C/C++ flags but add something only for this host:
# cc_extra		= ''

# Be picky about errors
# usepedantic		= 'no'

# Extra libraries
# sys_libs		= []

# Python libraries
# python_path		= '/sw/apps/python/x86_64/gcc-4.1.2/python-2.4.4/include/python2.4'
# python_lib_path	= '/sw/apps/python/x86_64/gcc-4.1.2/python-2.4.4/lib'
# python_libs		= 'python2.4'
# python_cmd		= 'python'

# Boost libraries
# boost_path		= '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.4.4/boost_1_33/include/boost-1_33'
# boost_lib_path	= '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.4.4/boost_1_33/lib'
# boost_libs		= ['boost_python-gcc']

# Specify whether or not to use VTK
# usevtk		= 'yes'

# NetCDF
# usenetcdf		= 'yes'
# netCDF_path		= '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/include'
# netCDF_lib_path	= '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/lib'
# netCDF_libs		= ['netcdf_c++', 'netcdf']

# MKL
# usemkl		= 'yes'
# mkl_path		= '/sw/sdev/cmkl/10.0.2.18/include'
# mkl_lib_path		= '/sw/sdev/cmkl/10.0.2.18/lib/em64t'
# mkl_libs		= ['mkl_solver', 'mkl_em64t', 'mkl_core', 'guide', 'pthread']

# UMFPACK (requires AMD and BLAS)
# useumfpack		= 'yes'
# ufc_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/include'
# umf_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/include'
# umf_lib_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/lib'
# umf_libs		= ['umfpack']
# amd_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/include'
# amd_lib_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/lib'
# amd_libs		= ['amd']
# blas_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/include'
# blas_lib_path		= '/sw/libs/umfpack/x86_64/gcc-4.1.2/umfpack-5.2/lib'
# blas_libs		= ['blas']

# OpenMP
# useopenmp		= 'yes'

# MPICH2 (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# mpi_path		= '/home/Work/InstallArea/mpich2-1.0.7/include'
# mpi_lib_path		= '/home/Work/InstallArea/mpich2-1.0.7/lib'
# mpi_libs		= ['mpich', 'rt']
# mpi_run		= 'mpirun -np 1'

# MPICH2 for jumpshot (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# mpi_path		= '/home/Work/InstallArea/mpich2-1.0.7/include'
# mpi_lib_path		= '/home/Work/InstallArea/mpich2-1.0.7/lib'
# mpi_libs		= ['lmpe', 'mpe', 'mpich', 'rt']
# mpi_run		= 'mpirun -np 1'

# ParMETIS (for use with MPI)
# useparmetis		= 'yes'
# parmetis_path		= '/home/Work/InstallArea/parmetis-3.1/include'
# parmetis_lib_path	= '/home/Work/InstallArea/parmetis-3.1/lib'
# parmetis_libs		= ['parmetis', 'metis']

# PAPI
# usepapi		= 'no'
# papi_path		= '/usr/include'
# papi_lib_path		= '/usr/lib'
# papi_libs		= ['papi']
# papi_instrument_solver	= 'no'

