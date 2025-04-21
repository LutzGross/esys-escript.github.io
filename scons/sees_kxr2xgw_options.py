
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
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
# sudo port install gmsh
# sudo port select --set mpi mpich-mp-fortran

# sudo port clang-14
# sudo port select --set pygments py310-pygments
# sudo rm -r  /opt/local/Library/Frameworks/Python.framework/Versions/3.10/lib/python3.10/site-packages/setuptools-63.2.0.dist-info

# sudo port install silo


# sudo port install lapack - build failed 28/12/22
# sudo port install py310-mpi4py

# "-I/opt/local/include/libomp -L/opt/local/lib/libomp -fopenmp" ???

from templates.homebrew_options import *
build_trilinos="True"

#cxx_extra = '-Wimplicit-function-declaration -Wno-string-concatenation -fopenmp' # --target arm64-apple-macosx13.0.0'
# ld_extra='-v'# -L/Users/uqlgross/PycharmProjects/esys-escript.github.io/build/darwin/escriptcore/src'


