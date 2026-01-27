
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# =============================================================================
# WINDOWS BUILD OPTIONS - DEPRECATED
# =============================================================================
#
# RECOMMENDED APPROACH: Use WSL2 (Windows Subsystem for Linux)
# -----------------------------------------------------------------------------
# Native Windows builds are complex and not actively tested. The recommended
# approach for Windows users is to use WSL2 with Ubuntu or Debian:
#
#   1. Install WSL2: https://docs.microsoft.com/en-us/windows/wsl/install
#   2. Install Ubuntu from the Microsoft Store
#   3. Inside WSL2, follow the standard Linux installation instructions
#   4. Use debian_options.py or ubuntu_options.py as your options file
#
# This provides a fully supported Linux environment on Windows with much
# simpler dependency management.
#
# =============================================================================
# NATIVE WINDOWS BUILD (unsupported, for reference only)
# =============================================================================
# The settings below are provided for reference only. Native Windows builds
# require manual configuration of all dependencies and are not tested.
#
# Copy this file to <hostname>_options.py, where <hostname> is your machine's
# short hostname, then customize to your needs.
#
# PREFIXES:
# There are two ways to specify where to find dependent headers and libraries:
# 1) If your installation follows the general scheme where headers are located
#    in <prefix>/include, and libraries in <prefix>/lib then it is sufficient
#    to specify this prefix, e.g. boost_prefix='C:/python'
# 2) Otherwise provide a list with two elements, where the first one is the
#    include path, and the second the library path, e.g.
#    boost_prefix=['C:/boost/include/boost1_84', 'C:/boost/lib']

# The options file version. SCons will refuse to build if there have been
# changes to the set of variables and your file has not been updated.
escript_opts_version = 203

# Installation prefix. Files will be installed in subdirectories underneath.
#prefix = 'C:/escript'

# C++ compiler command name or full path.
#cxx = 'cl'

# Additional flags to add to the C++ compiler
#cxx_extra = ''

# Additional flags to add to the linker
#ld_extra = ''

# Whether to treat compiler warnings as errors
#werror = False

# Whether to build a debug version
#debug = False

# Set to True to add flags that enable OpenMP parallelization
#openmp = True

# Flavour of MPI implementation
# Recognized values: 'none', 'MPICH', 'OPENMPI', 'INTELMPI', 'MPT'
#mpi = 'none'

# Prefix or paths to MPI headers and libraries
#mpi_prefix = 'C:/Program Files/Microsoft MPI'

# Prefix or paths to boost-python headers and libraries
#boost_prefix = 'C:/boost'

# boost-python library/libraries to link against
#boost_libs = ['boost_python312']

# Whether to use the netCDF library for dump file support
#netcdf = False

# Prefix or paths to netCDF headers and libraries
#netcdf_prefix = 'C:/netcdf'

# Whether to use UMFPACK
#umfpack = False

# Prefix or paths to UMFPACK headers and libraries
#umfpack_prefix = 'C:/SuiteSparse'

# Flavour of LAPACK implementation
# Recognized values: 'none', 'clapack', 'mkl'
#lapack = 'none'

# Whether to use LLNL's SILO library for output file support
#silo = False

# List of domain families to build
#domains = 'finley,ripley,speckley'

# Extra libraries to link with
#sys_libs = ['ws2_32']

# Tools for SCons
#tools_names = ['msvc']
