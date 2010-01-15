
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


# Configuration for guineapig (64-bit Intel running Debian)

# If you cannot use the default compiler flags set in SConstruct, then change them here
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
# cc_flags		= '-ansi'
# cc_optim		= '-O2'
# cc_debug		= '-g'
omp_optim		= '-fopenmp'
omp_debug		= '-fopenmp'
omp_libs		= []

# Use the default C/C++ flags but add something only for this host:
#cc_extra		= '-Wall -pedantic -isystem /usr/include/boost/ -isystem /usr/include/python2.5/ -Wno-sign-compare -Wno-long-long'
#cc_extra		= '-fopenmp'
ld_extra		= '-fopenmp'
cc_extra = '-isystem /usr/local/py2.6.2/silo4.7.2/include'

# Be picky about errors
# usepedantic		= 'no'

# Extra libraries
# sys_libs		= ['guide', 'pthread', 'stdc++']

# Python libraries
python_path		= '/usr/local/python2.6.2/include/python2.6'
python_lib_path		= '/usr/local/python2.6.2/lib'
python_libs		= 'python2.6'
# python_cmd		= 'python'

# Boost libraries
boost_path		= '/usr/local/py2.6.2/boost1.39.0/include/boost-1_39'
boost_lib_path		= '/usr/local/py2.6.2/boost1.39.0/lib'
boost_libs		= ['libboost_python-gcc44-mt']

# Specify whether or not to use VTK
# usevtk		= 'yes'

# NetCDF
#usenetcdf		= 'yes'
#netCDF_path		= '/usr/local/py2.6.2/netcdf4.0/include'
#netCDF_lib_path		= '/usr/local/py2.6.2/netcdf4.0/lib'
#netCDF_libs		= ['netcdf_c++', 'netcdf']

# Silo
usesilo         = 'yes'
#silo_path       = '/usr/local/py2.6.2/silo4.7.2/include'
silo_path       = '/usr/local/py2.6.2/silo4.7.2/include'
silo_lib_path   = '/usr/local/py2.6.2/silo4.7.2/lib'
silo_libs       = ['siloh5', 'hdf5']

# MKL
# usemkl		= 'yes'
# mkl_path		= '/sw/sdev/cmkl/10.0.2.18/include'
# mkl_lib_path		= '/sw/sdev/cmkl/10.0.2.18/lib/em64t'
# mkl_libs		= ['mkl_solver', 'mkl_em64t', 'mkl_core', 'guide', 'pthread']

# UMFPACK (requires AMD and BLAS)
useumfpack		= 'yes'
ufc_path		= '/usr/include/suitesparse'
umf_path		= '/usr/include/suitesparse'
umf_lib_path		= '/usr/lib'
umf_libs		= ['umfpack']
amd_path		= '/usr/include/suitesparse'
amd_lib_path		= '/usr/lib'
amd_libs		= ['amd']
blas_path		= '/usr/include'
blas_lib_path		= '/usr/lib'
blas_libs		= ['blas']

# OpenMP
useopenmp		= 'yes'

usempi 			= 'yes'
mpi_flavour 		= 'OPENMPI'
mpi_path		= '/usr/include/openmpi/'
mpi_lib_path		= '/usr/lib/openmpi/lib/'
mpi_libs		= ['libmpi','libmpi_cxx']

# MPICH2 (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# usempi		= 'no'
# mpi_path		= '/home/Work/InstallArea/mpich2-1.0.7/include'
# mpi_lib_path		= '/home/Work/InstallArea/mpich2-1.0.7/lib'
# mpi_libs		= ['mpich', 'rt']
# mpi_flavour		= "MPICH"

# MPICH2 for jumpshot (to run Escript use: module load mpich2/gcc-4.1.2/mpich2-1.0.7)
# mpi_path		= '/home/Work/InstallArea/mpich2-1.0.7/include'
# mpi_lib_path		= '/home/Work/InstallArea/mpich2-1.0.7/lib'
# mpi_libs		= ['lmpe', 'mpe', 'mpich', 'rt']
# mpi_flavour		= "MPICH"

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

