
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
#   brew install hdf5 suite-sparse
#   brew install open-mpi  # Optional, for MPI support
#   pip3 install numpy scipy matplotlib
#
# For LLVM compiler support (recommended for OpenMP):
#   echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc
#   export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
#   export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"

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

# MPI configuration - disabled by default on macOS
# To enable MPI, install openmpi: brew install open-mpi
# Then set mpi = 'OPENMPI' and configure mpi_prefix
mpi = 'none'

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

# Domain configuration - oxley excluded (requires p4est)
domains = ('finley', 'ripley', 'speckley')

# OpenMP configuration
openmp = True
omp_flags = ["-fopenmp"]
omp_ldflags = ["-fopenmp"]

# Boost configuration
boost_prefix = HOMEBREW_PREFIX
boost_libs = [f'boost_python3{subversion}']
compression_libs = ['boost_iostreams']

# HDF5 configuration
# Install with: brew install hdf5
hdf5 = True
hdf5_prefix = HOMEBREW_PREFIX
hdf5_libs = ['hdf5_cpp', 'hdf5']

# SILO configuration (for VisIt visualization)
# Install with: brew install silo
silo = True
silo_prefix = HOMEBREW_PREFIX
silo_libs = ['siloh5', 'hdf5']

# NetCDF configuration
# Install with: brew install netcdf netcdf-cxx4
netcdf = True
netcdf_prefix = HOMEBREW_PREFIX
netcdf_libs = ['netcdf-cxx4', 'netcdf']

# UMFPACK direct solver
# Install with: brew install suite-sparse
umfpack = True
umfpack_prefix = [HOMEBREW_PREFIX + "/include/suitesparse", HOMEBREW_PREFIX + "/lib/"]

# LAPACK configuration - auto-detect
lapack = 'auto'

# Trilinos - build from bundled source
build_trilinos = 'make'

# zlib - typically available on macOS
zlib = True
zlib_prefix = HOMEBREW_PREFIX
zlib_libs = ['z']

# Linker flags
ld_extra = ["-L/opt/homebrew/opt/llvm/lib", "-lz", "-Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++"]

# Clean up temporary variables
del HOMEBREW_PREFIX, p, subversion, revversion, t
