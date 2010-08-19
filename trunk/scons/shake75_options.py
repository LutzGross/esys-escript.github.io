
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


# Configuration for shake75 (32-bit Intel Core 2 running Debian Linux)
# If you cannot use the default compiler flags set in SConstruct,
# then change them here
#cc="gcc-4.4"
#cxx="g++-4.4"
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
# cc_flags      = '-ansi'
# cc_optim      = '-O2'
# usedebug      = 'no'
# cc_debug      = '-g'
omp_optim       = '-fopenmp'
omp_debug       = ''
omp_libs        = []

# Use the default C/C++ flags but add something only for this host:
cc_extra        = '-Wall -mmmx -msse'
cxx_extra       = '-Wall -mmmx -msse'
ld_extra        = '-fopenmp'

# Be picky about errors
# usepedantic   = 'no'

# Extra libraries
# sys_libs      = []

# Python libraries
# python_path       = '/usr/lib/python2.6'
# python_lib_path   = '/usr/lib'
# python_libs       = 'python2.6'
# python_cmd        = 'python'

# Boost libraries
boost_path          = '/usr/include/'
# boost_lib_path    = '/usr/lib'
boost_libs          = ['boost_python']

# Specify whether or not to use VTK
# usevtk            = 'yes'

# NetCDF
#usenetcdf          = 'yes'
#netCDF_path        = '/usr/include/'
#netCDF_lib_path    = '/usr/lib/'
#netCDF_libs        = ['netcdf_c++', 'netcdf']

# MKL
# usemkl        = 'yes'
# mkl_path      = '/usr/include'
# mkl_lib_path  = '/usr/lib'
# mkl_libs      = ['mkl_solver', 'mkl_em64t', 'mkl_core', 'guide', 'pthread']

# UMFPACK (requires AMD and BLAS)
useumfpack    = 'yes'
# ufc_path      = '/usr/include'
# umf_path      = '/usr/include'
# umf_lib_path  = '/usr/lib'
# umf_libs      = ['umfpack']
# amd_path      = '/usr/include'
# amd_lib_path  = '/usr/lib'
# amd_libs      = ['amd']
# blas_path     = '/usr/include'
# blas_lib_path = '/usr/lib'
# blas_libs     = ['blas']

# Silo
usesilo         = 'yes'
silo_path       = '/usr/local/include/'
silo_lib_path   = '/usr/local/lib/'
silo_libs       = ['siloh5', 'hdf5']

# VisIt
usevisit        = 'yes'
visit_path      = '/opt/visit/2.1.0/linux-intel/libsim/V2/include'
visit_lib_path  = '/opt/visit/2.1.0/linux-intel/libsim/V2/lib'
#visit_path      = '/home/caltinay/src/visit-trunk/src/include/visit'
#visit_lib_path  = '/home/caltinay/src/visit-trunk/src/lib'

# OpenMP
useopenmp       = 'yes'

# MPI
usempi          = 'yes'
mpi_path        = '/usr/include/mpi'
mpi_lib_path    = '/usr/lib/openmpi/lib'
mpi_libs        = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
mpi_flavour     = 'OPENMPI'

# ParMETIS (for use with MPI)
#useparmetis        = 'yes'
#parmetis_path      = '/usr/local/include'
#parmetis_lib_path  = '/usr/local/lib'
#parmetis_libs      = ['parmetis', 'metis']

# PAPI
# usepapi                = 'no'
# papi_path              = '/usr/include'
# papi_lib_path          = '/usr/lib'
# papi_libs              = ['papi']
# papi_instrument_solver = 'no'

