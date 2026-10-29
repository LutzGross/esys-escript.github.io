
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

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
#cxx_extra = '-Wno-literal-suffix'
openmp = True
mpi='OPENMPI'

prefix = '/app'

cxx_extra = '-w -O3 -march=native'

#pythoncmd='/app/bin/python3'
pythonlibpath = '/app/lib'
pythonincpath = '/usr/include/python3.14'

boost_prefix = ['/app/include','/app/lib']
boost_libs = ['boost_python314','boost_iostreams','boost_random']

domains = ['finley','ripley','speckley']

hdf5_prefix = ['/app/include','/app/lib']
hdf5_libs = ['hdf5_cpp','hdf5','hdf5_hl_cpp']

mpi_prefix = ['/app/include','/app/lib']
mpi_libs = ['mpi']
mpi4py = True

lapack=1
lapack_prefix = ['/app/include/', '/app/lib64/']
lapack_libs = ['lapacke','lapack', 'cblas', 'blas']

ld_extra=''
paso=1
p4est=0

netcdf = True
netcdf_prefix = ['/app/include', '/app/lib/']
netcdf_libs = ['netcdf_c++4', 'netcdf']

zlib = True
zlib_libs = ['z']

trilinos = True
trilinos_prefix = ['/app/include/', '/app/lib']
trilinos_make_sh = 'tools/flatpak/flatpak_trilinos.sh'

umfpack = False
umfpack_prefix = ['/app/include','/app/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

silo = True
silo_prefix = ['/app/include', '/app/lib']
silo_libs = ['siloh5', 'hdf5']

sympy = True

visit = 0
werror = 0


# boost-python library/libraries to link against
