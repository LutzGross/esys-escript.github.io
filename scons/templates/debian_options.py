##############################################################################
#
# 
#
##############################################################################

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
sysheaderopt = False
openmp = True
paso=True

cxx_extra= ['-O3', '-fdiagnostics-color=always', '-fstack-protector-strong',  '-Wformat', '-Werror=format-security' ]


pythoncmd = 'python3'
import subprocess
p=subprocess.run([pythoncmd,'-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1]
print("python subversion = ",subversion)
pythonlibname = 'python3.%s'%subversion
pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath =  '/usr/include/python3.%s'%subversion

boost_libs=['boost_python3%s'%subversion,'boost_numpy3%s'%subversion,'boost_random']
boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
werror=0

mpi='NO'
mpi_prefix = ['/usr/include/x86_64-linux-gnu/openmpi', '/usr/lib/x86_64-linux-gnu/openmpi' ]
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']

# METIS configuration - graph partitioning for Trilinos
# Install with: sudo apt-get install libmetis-dev
metis = True
metis_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
metis_libs = ['metis']

# ParMETIS configuration - parallel graph partitioning (requires MPI)
# Install with: sudo apt-get install libparmetis-dev
parmetis = True
parmetis_prefix = ['/usr/include/parmetis', '/usr/lib/x86_64-linux-gnu']
parmetis_libs = ['parmetis']

umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
#umfpack_libs = ['umfpack', 'blas', 'amd']

silo = True
silo_libs = ['siloh5', 'hdf5']
silo_prefix=[ '/usr/include' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']

hdf5 = True
hdf5_libs = ['hdf5_serial_cpp']
hdf5_prefix=[ '/usr/include/hdf5/serial' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']

# LAPACK configuration - uses LAPACKE (modern C interface)
# Install with: sudo apt-get install liblapacke-dev
lapack = 'auto'  # Auto-detect LAPACKE
lapack_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
lapack_libs = ['lapacke']

# zlib configuration - required by p4est for oxley domain
# Install with: sudo apt-get install zlib1g-dev
zlib = True
zlib_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
zlib_libs = ['z']
