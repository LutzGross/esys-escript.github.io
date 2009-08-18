
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


# Configuration for Cognac (SGI Altix)

#  Assumed modules:
#
# module load intel-cc/10.0.023
# module load python/2.4.4/icc10.0.023
# module load numarray/1.5.2/python2.4.4_icc10.0.023/icc10.0.023
# module load Mesa/7.0.2
# module load VTK/5.0.3_mesa
# module load boost/1.33.1/python2.4.4_icc10.0.023/numarray1.5.2_icc10.0.023/icc10.0.023
# module load netcdf/3.6.2
# module load intel-mkl/9.1.018
# module load scsl/1.6.1.0
# module load scons/0.97
# module load subversion/1.4.6
#

python_version="2.4"
python_installation="2.4.4/icc10.0.023"
boost_version="1_33_1"
boost_installation="1.33.1/python2.4.4_icc10.0.023/numarray1.5.2_icc10.0.023/icc10.0.023"

# If you cannot use the default compiler flags set in SConstruct, then change them here
# C/C++ Compiler flags (always use cc_flags and either cc_optim or cc_debug)
cc_flags		= '-fPIC -ansi -wd161 -w1 -DBLOCKTIMER -DCORE_ID1'
# cc_optim		= '-O2'
# cc_debug		= '-g'
# omp_optim		= '-openmp'
# omp_debug		= '-openmp'
# omp_libs		= ['guide']

# Use the default C/C++ flags but add something only for this host:
# cc_extra		= ''
# ld_extra		= ''

# Be picky about errors
# usepedantic		= 'no'

# Extra libraries
# sys_libs		= []

# Python libraries
python_path		= '/opt/python/'+python_installation+'/include/python'+python_version
python_lib_path		= '/opt/python/'+python_installation+'/lib'
python_libs		= ['python'+python_version]
# python_cmd		= 'python'

# Boost libraries
boost_path		= '/opt/boost/'+boost_installation+'/include/boost-'+boost_version
boost_lib_path		= '/opt/boost/'+boost_installation+'/lib'
boost_libs		= ['boost_python-il-mt-1_33_1']

# Specify whether or not to use VTK
# usevtk		= 'yes'

# NetCDF
usenetcdf		= 'yes'
netCDF_path		= '/opt/netcdf/3.6.2/include'
netCDF_lib_path		= '/opt/netcdf/3.6.2/lib'
# netCDF_libs		= ['netcdf_c++', 'netcdf']

# MKL
# usemkl		= 'no'
# mkl_path		= '/opt/intel_mkl/9.1.018/include'
# mkl_lib_path		= '/opt/intel_mkl/9.1.018/lib/64'
# mkl_libs		= ['mkl_solver', 'mkl_lapack', 'mkl_ipf']

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
usempi			= 'no'
mpi_path		= '/usr/include'
mpi_lib_path		= '/usr/lib'
mpi_libs		= ['mpi']
mpi_flavour		= 'MPT'

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

