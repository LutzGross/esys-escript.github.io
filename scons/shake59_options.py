
#          Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#

# Configuration for shake59 (Windows)

# If you cannot use the default compiler flags set in SConstruct, then change them here
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
# cc_flags		= '-ansi'
# cc_optim		= '-O2'
# cc_debug		= '-g'
# omp_optim		= '-openmp'
# omp_debug		= '-openmp'
# omp_libs		= ['guide']

# Use the default C/C++ flags but add something only for this host:
# cc_extra		= ''

# Be picky about errors
# usepedantic		= 'no'

# Extra libraries
# sys_libs		= ['guide', 'pthread', 'stdc++']

# Python libraries
python_path		= 'C:/python25/include'
python_lib_path		= 'C:/python25/libs'
python_libs		= ['python25']
# python_cmd		= 'python'

# Boost libraries
boost_path		= 'C:/Documents and Settings/ksteube/Desktop/Temp/boost_1_33'
boost_lib_path		= 'C:/Documents and Settings/ksteube/Desktop/Temp/boost_1_33/windows_binary/lib'
boost_libs		= ['boost_python-vc71-mt-1_33']

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

# MPI MPT (no module load required)
# usempi		= 'no'
# mpi_path		= '/usr/include'
# mpi_lib_path		= '/usr/lib64'
# mpi_libs		= ['mpi']
# mpi_run		= 'mpirun -np 1'

# MPICH2 (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# mpi_path		= '/sw/libs/mpich2/x86_64/gcc-4.1.2/mpich2-1.0.7/include'
# mpi_lib_path		= '/sw/libs/mpich2/x86_64/gcc-4.1.2/mpich2-1.0.7/lib'
# mpi_libs		= ['mpich', 'rt']
# mpi_run		= 'mpirun -np 1'

# MPICH2 for jumpshot (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# mpi_path		= '/sw/libs/mpich2/x86_64/gcc-4.1.2/mpich2-1.0.7/include'
# mpi_lib_path		= '/sw/libs/mpich2/x86_64/gcc-4.1.2/mpich2-1.0.7/lib'
# mpi_libs		= ['lmpe', 'mpe', 'mpich', 'rt']
# mpi_run		= 'mpirun -np 1'

# Intel MPI (to run Escript use: module load intel-mpi/3.1.038)
# mpi_path		= '/sw/sdev/intel/mpi/3.1.038/x86_64/include64'
# mpi_lib_path		= '/sw/sdev/intel/mpi/3.1.038/x86_64/lib64'
# mpi_libs		= ['mpi']
# mpi_run		= 'mpirun -np 1'

# mvapich (to run Escript use: module load mvapich/mvapich-1.0.1)
# mpi_path		= '/usr/diags/mpi/mvapich/intel/include'
# mpi_lib_path		= '/usr/diags/mpi/mvapich/intel/lib/shared'
# mpi_libs		= ['mpich']
# mpi_run		= 'mpirun -np 1'

# OpenMPI (to run Escript use: module load openmpi/gcc-4.1.2/openmpi-1.2.6) (This doesn't compile)
# mpi_path		= '/sw/libs/openmpi/x86_64/gcc-4.1.2/openmpi-1.2.6/include'
# mpi_lib_path		= '/sw/libs/openmpi/x86_64/gcc-4.1.2/openmpi-1.2.6/lib'
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

