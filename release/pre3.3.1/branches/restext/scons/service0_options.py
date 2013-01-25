
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


# Configuration for Savanna (SGI ICE 8200 running SUSE Linux)

# Append environment variables which need to be passed through scons to 
# other tools
env_export=['INTEL_LICENSE_FILE']

# If you cannot use the default compiler flags set in SConstruct, then change them here
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
# cc_flags		= '-ansi'
# cc_optim		= '-O2'
# cc_debug		= '-g'
omp_optim		= '-openmp  -openmp-report2'
# omp_debug		= '-openmp'
# omp_libs		= ['guide']

# Use the default C/C++ flags but add something only for this host:
# cc_extra		= ''
ld_extra		= '-shared-intel'	# Fix warning about feupdate in icc v10

# Be picky about errors
# usepedantic		= 'no'

# Extra libraries
# sys_libs		= []

# Python libraries
python_path		= '/sw/apps/python/x86_64/gcc-4.1.2/python-2.6.2/include/python2.6'
python_lib_path		= '/sw/apps/python/x86_64/gcc-4.1.2/python-2.6.2/lib'
# python_libs		= ['python2.4']
# python_cmd		= 'python'

# Boost libraries
boost_path		= '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.6.2/boost_1_39_0/include/boost-1_39'
boost_lib_path		= '/sw/libs/boost/x86_64/gcc-4.1.2/python-2.6.2/boost_1_39_0/lib'
boost_libs		= ['boost_python-gcc41-mt']

# Specify whether or not to use VTK
# usevtk		= 'yes'

# NetCDF
usenetcdf		= 'yes'
netCDF_path		= '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/include'
netCDF_lib_path		= '/sw/libs/netcdf/x86_64/gcc-4.1.2/netcdf-3.6.2/lib'
# netCDF_libs		= ['netcdf_c++', 'netcdf']

# MKL
usemkl		= 'yes'
mkl_path		= '/sw/sdev/intel/cmkl/10.1.0.015/include'
mkl_lib_path		= '/sw/sdev/intel/cmkl/10.1.0.015/lib/em64t'
mkl_libs                =  [ "mkl_core", "mkl_intel_lp64",  "mkl_intel_thread", "mkl_lapack", 'guide', 'pthread' , " mkl_mc", "mkl_def"]
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

# ParMETIS (for use with MPI) (these setting my be overwritten depending on the mpi version used)
useparmetis		= 'yes'
parmetis_path		= '/sw/libs/parmetis/x86_64/gcc-4.1.2/intelmpi/parmetis-3.1/include'
parmetis_lib_path	= '/sw/libs/parmetis/x86_64/gcc-4.1.2/intelmpi/parmetis-3.1/lib'
parmetis_libs		= ['parmetis', 'metis']
# Silo
# usesilo		= 'yes'
silo_path		= '/sw/libs/silo/x86_64/gcc-4.1.2/silo-4.6.1/include'
silo_lib_path		= '/sw/libs/silo/x86_64/gcc-4.1.2/silo-4.6.1/lib'
# silo_libs		= ['siloh5', 'hdf5']


# OpenMP
useopenmp		= 'yes'

usempi		= 'yes'
# MPI MPT (no module load required)
# mpi_path		= '/usr/include'
# mpi_lib_path		= '/usr/lib64'
# mpi_libs		= ['mpi']
# mpi_flavour     = "MPT"

# MPICH2 (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# mpi_path		= '/sw/libs/mpich2/x86_64/gcc-4.1.2/mpich2-1.0.7/include'
# mpi_lib_path		= '/sw/libs/mpich2/x86_64/gcc-4.1.2/mpich2-1.0.7/lib'
# mpi_libs		= ['mpich', 'rt']
# mpi_flavour		= 'MPICH2'

# Intel MPI (to run Escript use: module load intel-mpi/3.2.0.011
mpi_path		= '/sw/sdev/intel/mpi/3.2.0.011/x86_64/include64'
mpi_lib_path		= '/sw/sdev/intel/mpi/3.2.0.011/x86_64/lib64'
mpi_libs		= ['mpi']
mpi_flavour		= 'INTELMPI'

# ParMETIS setting needs to be overwritten in intel-mpi is used
# parmetis_path = "/sw/libs/parmetis/x86_64/gcc-4.1.2/intelmpi/parmetis-3.1/include"
# parmetis_lib_path   = "/sw/libs/parmetis/x86_64/gcc-4.1.2/intelmpi/parmetis-3.1/lib"

# mvapich (to run Escript use: module load mvapich/mvapich-1.0.1)
# mpi_path		= '/usr/diags/mpi/mvapich/intel/include'
# mpi_lib_path		= '/usr/diags/mpi/mvapich/intel/lib/shared'
# mpi_libs		= ['mpich']
# mpi_flavour		= 'MPICH'

# OpenMPI (to run Escript use: module load openmpi/gcc-4.1.2/openmpi-1.2.6) (This doesn't compile)
# mpi_path		= '/sw/libs/openmpi/x86_64/gcc-4.1.2/openmpi-1.2.6/include'
# mpi_lib_path		= '/sw/libs/openmpi/x86_64/gcc-4.1.2/openmpi-1.2.6/lib'
# mpi_libs		= ['mpi']
# mpi_flavour		= 'OPENMPI'

# PAPI
# usepapi		= 'no'
# papi_path		= '/usr/include'
# papi_lib_path		= '/usr/lib64'
# papi_libs		= ['papi']
# papi_instrument_solver	= 'no'

