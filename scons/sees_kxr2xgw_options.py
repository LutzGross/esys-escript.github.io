
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
#
# Prerequisites (install via Homebrew):
#   brew install python3 scons cmake llvm
#   brew install boost boost-python3
#   brew install hdf5 suite-sparse netcdf netcdf-cxx4
#   brew install open-mpi mpi4py
#   pip3 install numpy scipy sympy matplotlib
#
# For LLVM compiler support (recommended for OpenMP):
#   echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc
#
# See scons/templates/homebrew_options.py for full configuration details

from templates.homebrew_options import *

# Override to always build Trilinos from source
build_trilinos = 'make'

# Skip library link checks on macOS (Python 3.14 compatibility)
# Libraries are verified via Homebrew and will link correctly during build
skip_link_checks = True
