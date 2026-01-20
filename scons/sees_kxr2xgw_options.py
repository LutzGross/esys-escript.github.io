
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

# This is a template configuration file for escript on Mac OS X using ports.
# Refer to README_FIRST for usage instructions.
#
# install Xcode + line commands (start Xcode to accept licence conditions)


# sudo port install clang-14
# sudo port select --set pygments py310-pygments
# sudo port install scons
# sudo port select --set python python310
# sudo port select --set python3 python
# sudo port install boost-numpy
# sudo port select --set cython cython310
# sudo port install netcdf-cxx4
# sudo port install SuiteSparse
# sudo port install py310-matplotlib-inline
# sudo port install py310-sympy
# port select --set py-sympy py310-sympy
# sudo port install py310-matplotlib
# sudo port install py310-scipy
# sudo port select --set pip pip310
# sudo port select --set pip2 pip310
# sudo port select --set mpi mpich-mp-fortran

# sudo port clang-14
# sudo port select --set pygments py310-pygments
# sudo rm -r  /opt/local/Library/Frameworks/Python.framework/Versions/3.10/lib/python3.10/site-packages/setuptools-63.2.0.dist-info

# sudo port install silo


# sudo port install lapack - build failed 28/12/22
# sudo port install py310-mpi4py

# "-I/opt/local/include/libomp -L/opt/local/lib/libomp -fopenmp" ???

from templates.homebrew_options import *
build_trilinos='make'

#cxx_extra = '-Wimplicit-function-declaration -Wno-string-concatenation -fopenmp' # --target arm64-apple-macosx13.0.0'
# ld_extra='-v'# -L/Users/uqlgross/PycharmProjects/esys-escript.github.io/build/darwin/escriptcore/src'


