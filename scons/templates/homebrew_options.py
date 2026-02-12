
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

# This is a template configuration file for escript on macOS with Homebrew.
# Refer to installation.md for usage instructions.
#
# Prerequisites:
#   brew install python3 scons cmake llvm
#   brew install boost boost-python3
#   brew install hdf5 suite-sparse netcdf netcdf-cxx4 zlib metis
#   brew install open-mpi mpi4py
#   brew install numpy scipy python-matplotlib
#
# Optional (sympy not available in Homebrew):
#   For symbolic math support, install sympy in a virtual environment:
#   python3 -m venv ~/venv && source ~/venv/bin/activate && pip install sympy
#
# For LLVM compiler support (recommended for OpenMP):
#   echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc
#   export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
#   export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
#
# For no-MPI build, use homebrew_nompi_options.py instead.

import os
import subprocess

escript_opts_version = 203

# Homebrew prefix (Apple Silicon: /opt/homebrew, Intel: /usr/local)
HOMEBREW_PREFIX = '/opt/homebrew'

# Python configuration - auto-detect version
pythoncmd = os.path.join(HOMEBREW_PREFIX, 'bin', 'python3')
p = subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
print(p.stdout)
subversion = p.stdout.split(' ')[1].split('.')[1].strip()
revversion = p.stdout.split(' ')[1].split('.')[2].strip()
print(f'Compiled for Python 3.{subversion}.{revversion}.')

# Get Python include and library paths
p = subprocess.run([pythoncmd + f".{subversion}-config", '--includes'], capture_output=True, text=True)
for t in p.stdout.split(" "):
    if t.startswith("-I"):
        pythonincpath = t[2:]
        break

p = subprocess.run([pythoncmd + f".{subversion}-config", '--ldflags'], capture_output=True, text=True)
for t in p.stdout.split(" "):
    if t.startswith("-L"):
        pythonlibpath = t[2:]
        break

# Set SDKROOT for macOS builds
try:
    p = subprocess.run(['xcrun', '--show-sdk-path'], capture_output=True, text=True)
    if p.stdout.strip():
        os.putenv("SDKROOT", p.stdout.strip())
except FileNotFoundError:
    print("xcrun not found. SDKROOT is not set.")

# Compiler configuration - use LLVM for OpenMP support
# Install with: brew install llvm
cxx = '/opt/homebrew/opt/llvm/bin/clang++'
cc = '/opt/homebrew/opt/llvm/bin/clang'

cxx_extra = ["-Wno-error=uninitialized",
             "-Wno-error=implicit-function-declaration",
             "-Wno-error=unused-but-set-variable",
             "-Wno-error=return-stack-address",
             "-Wno-error=inconsistent-missing-override",
             "-Wno-error=unused-function"]

cc_extra = ["-Wno-error=unused-but-set-variable",
            "-Wno-error=deprecated-non-prototype"]
cc_optim = ["-O3"]

print(f"C++ compiler is {cxx}")

# OpenMP configuration
openmp = True
omp_flags = ["-fopenmp"]
omp_ldflags = ["-fopenmp"]

# MPI configuration
# Install with: brew install open-mpi mpi4py
mpi = 'OPENMPI'
# OpenMPI 5.x uses the standard include/lib directories (not openmpi subdirectories)
mpi_prefix = [HOMEBREW_PREFIX + '/include', HOMEBREW_PREFIX + '/lib']
# OpenMPI 5.x dropped libmpi_cxx - C++ bindings are now in libmpi
mpi_libs = ['mpi']
mpi4py = True

# Solver configuration
paso = True
build_trilinos = 'check'

# Domain configuration - oxley excluded (requires p4est)
domains = ('finley', 'ripley', 'speckley')

# Boost configuration
boost_prefix = HOMEBREW_PREFIX
boost_libs = [f'boost_python3{subversion}']
compression_libs = ['boost_iostreams']

# HDF5 configuration
hdf5 = True
hdf5_prefix = HOMEBREW_PREFIX
hdf5_libs = ['hdf5_cpp', 'hdf5']

# SILO configuration (for VisIt visualization)
silo = True
silo_prefix = HOMEBREW_PREFIX
silo_libs = ['siloh5', 'hdf5']

# NetCDF configuration
netcdf = True
netcdf_prefix = HOMEBREW_PREFIX
netcdf_libs = ['netcdf-cxx4', 'netcdf']

# UMFPACK direct solver
umfpack = True
umfpack_prefix = [HOMEBREW_PREFIX + "/include/suitesparse", HOMEBREW_PREFIX + "/lib/"]

# LAPACK configuration - auto-detect
lapack = 'auto'

# zlib - keg-only in Homebrew (macOS provides system version)
# Use Homebrew version for consistency: brew install zlib
zlib = True
zlib_prefix = [HOMEBREW_PREFIX + '/opt/zlib/include', HOMEBREW_PREFIX + '/opt/zlib/lib']
zlib_libs = ['z']

# METIS graph partitioning library (used by Trilinos)
# Install with: brew install metis
metis = True
metis_prefix = [HOMEBREW_PREFIX + '/include', HOMEBREW_PREFIX + '/lib']
metis_libs = ['metis']

# ParMETIS - not available in Homebrew, would need to build from source
parmetis = False

# Linker flags - include LLVM libc++ for macOS arm64
ld_extra = ["-L/opt/homebrew/opt/llvm/lib", "-L/opt/homebrew/opt/llvm/lib/c++", "-lz", "-Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++"]

# Optional features
# sympy not available in Homebrew - requires virtual environment installation
sympy = False

# Clean up temporary variables
del HOMEBREW_PREFIX, p, subversion, revversion, t
