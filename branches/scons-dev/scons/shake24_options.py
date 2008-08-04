
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
# omp_libs		= ['guide']

# Be picky about errors
usepedantic		= 'no'

# Extra libraries
# sys_libs		= []

# Python libraries
python_path		= '/usr/include/python2.5'
python_lib_path		= '/usr/lib'
# python_libs		= ['python2.5']
# python_cmd		= 'python'

# Boost libraries
boost_path		= '/opt/local/include/boost-1_35'
boost_lib_path		= '/opt/local/lib'
boost_libs		= ['boost_python']

# Specify whether or not to use VTK
# usevtk		= 'yes'

# NetCDF
# usenetcdf		= 'yes'
netCDF_path		= '/opt/local/include'
netCDF_lib_path		= '/opt/local/lib'
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

# MPI MPT (no module load required)
# usempi		= 'no'
# mpi_path		= '/usr/include'
# mpi_lib_path		= '/usr/lib64'
# mpi_libs		= ['mpi']
# mpi_run		= 'mpirun -np 1'

# ParMETIS (for use with MPI)
# useparmetis		= 'yes'
# parmetis_path		= '/sw/libs/parmetis/x86_64/gcc-4.1.2/parmetis-3.1/include'
# parmetis_lib_path	= '/sw/libs/parmetis/x86_64/gcc-4.1.2/parmetis-3.1/lib'
# parmetis_libs		= ['parmetis', 'metis']

# PAPI
# usepapi		= 'no'
# papi_path		= '/usr/include'
# papi_lib_path		= '/usr/lib64'
# papi_libs		= ['papi']
# papi_instrument_solver	= 'no'

