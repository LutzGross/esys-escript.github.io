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

pythoncmd = 'python3'
import subprocess
p=subprocess.run([pythoncmd,'-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1]
print("python subversion = ",subversion)
#pythonlibname = 'python3.%s'%subversion
#pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = '/usr/include/python3.%s'%subversion

boost_libs=['boost_python3%s'%subversion,'boost_numpy3%s'%subversion,'boost_random']
#boost_libs='boost_python3%s'%subversion
#boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
werror=0

mpi = 'OPENMPI'
mpi_prefix = '/usr/lib/openmpi'
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
parmetis = True


umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
#umfpack_libs = ['umfpack', 'blas', 'amd']

silo = True
#silo_libs = ['siloh5', 'hdf5']
silo_prefix=[ '/usr/include' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']

hdf5 = True
hdf5_libs = ['hdf5_serial_cpp']
hdf5_prefix=[ '/usr/include/hdf5/serial' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']
