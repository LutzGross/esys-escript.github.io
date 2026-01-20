
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

# This is a template configuration file for escript on OS X homebrew.
# Refer to README_FIRST for usage instructions.


#If you need to have llvm first in your PATH, run:
#  echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc
#
#For compilers to find llvm you may need to set:
#  export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
#  export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"

import os
escript_opts_version = 203

HOMEBREW_PREFIX = '/opt/homebrew'

pythoncmd = os.path.join(HOMEBREW_PREFIX, 'bin', 'python3')
import subprocess
p=subprocess.run([pythoncmd, '-V'], capture_output=True, text=True)
print(p.stdout)
subversion=p.stdout.split(' ')[1].split('.')[1].strip()
revversion=p.stdout.split(' ')[1].split('.')[2].strip()
print(f'Compiled for Python 3.{subversion}.{revversion}.')

p=subprocess.run([pythoncmd+f".{subversion}-config", '--includes'], capture_output=True, text=True)
for t in p.stdout.split(" ") :
    if t.startswith("-I"):
        pythonincpath=p.stdout.split(" ")[0][2:]

p=subprocess.run([pythoncmd+f".{subversion}-config", '--ldflags'], capture_output=True, text=True)
for t in p.stdout.split(" "):
    if t.startswith("-L") :
        pythonlibpath=p.stdout.split(" ")[0][2:]

try:
    p=subprocess.run(['xcrun', '--show-sdk-path'], capture_output=True, text=True)
    if p.stdout.strip():
        os.putenv("SDKROOT", p.stdout.strip())
except FileNotFoundError:
    print("xcrun not found. SDKROOT is not set.")
#pythonpath = os.path.join(HOMEBREW_PREFIX,  "Cellar", f"python@3.{subversion}",
#                             f"3.{subversion}.{revversion}", "Frameworks", "Python.framework",
#                             "Versions", f"3.{subversion}" )
#
#pythonlibpath = os.path.join(pythonpath,  "lib" )
#                             #f"python3.{subversion}", f"config-3.{subversion}-darwin")
#pythonincpath = os.path.join(pythonpath, "include", f"python3.{subversion}" )

mpi = 'none'

cxx = '/opt/homebrew/opt/llvm/bin/clang++'
cxx_extra = ["-Wno-error=uninitialized",
             "-Wno-error=implicit-function-declaration",
             "-Wno-error=unused-but-set-variable",
             "-Wno-error=return-stack-address",
             "-Wno-error=inconsistent-missing-override",
             "-Wno-error=unused-function"]

cc = '/opt/homebrew/opt/llvm/bin/clang'
cc_extra = ["-Wno-error=unused-but-set-variable",
            "-Wno-error=deprecated-non-prototype" ]
cc_optim     = ["-O3" ]

print(f"c++ compiler is {cxx}")
domains = ('finley', 'ripley','speckley')

cppunit_prefix = HOMEBREW_PREFIX

openmp = True
omp_flags = ["-fopenmp"]
omp_ldflags = ["-fopenmp"]

boost_prefix = HOMEBREW_PREFIX
boost_libs = [f'boost_python3{subversion}']
compression_libs = ['boost_iostreams']

#boost_libs = [f'boost_python3{subversion}-mt']
#compression_libs = ['boost_iostreams-mt']
#============================================================
#netcdf = 4
#netcdf_prefix = HOMEBREW_PREFIX
#netcdf_libs=['netcdf-cxx4', 'netcdf']
hdf5 = True
hdf5_prefix = HOMEBREW_PREFIX
hdf5_libs=['hdf5_cpp', 'hdf5']
#===========================================================
silo = True
silo_prefix = HOMEBREW_PREFIX
silo_libs = ['siloh5', 'hdf5' ]

zlib = False
zlib_prefix = HOMEBREW_PREFIX+"/Cellar/zlib/1.3"
zlib_libs = [ 'zlib']

lapack =False
lapack_prefix = os.path.join(HOMEBREW_PREFIX, 'Cellar', 'lapack', '3.12.0' )

umfpack = True
umfpack_prefix = [ HOMEBREW_PREFIX+"/include/suitesparse", HOMEBREW_PREFIX+"/lib/" ]

build_trilinos = 'make'

ld_extra = ["-L/opt/homebrew/opt/llvm/lib", "-lz", "-Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++" ]


# clean up:
del HOMEBREW_PREFIX, p, subversion, revversion, t
