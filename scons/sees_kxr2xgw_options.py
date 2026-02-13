
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

# Configuration file for escript on macOS with Homebrew
# Host: sees-kxr2xgw (Apple Silicon)
# Build: NO-MPI (single-process only)
#
# Prerequisites (install via Homebrew):
#   brew install python3 scons cmake llvm
#   brew install boost boost-python3
#   brew install hdf5 suite-sparse netcdf netcdf-cxx4
#   pip3 install numpy scipy sympy matplotlib
#
# For LLVM compiler support (recommended for OpenMP):
#   echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc
#
# See scons/templates/homebrew_nompi_options.py for full configuration details

from templates.homebrew_nompi_options import *

# Build Trilinos without MPI (serial version)
trilinos = True
build_trilinos = 'make'
