##############################################################################
#
# 
#
##############################################################################

# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

escript_opts_version = 203
openmp = True
paso=True

pythoncmd = 'python3'
import subprocess
p=subprocess.run([pythoncmd,'-V'], capture_output=True, text=True)
subversion=p.stdout.split(' ')[1].split('.')[1]
print("subversion = ",subversion)
pythonlibname = 'python3.%s'%subversion
pythonlibpath = '/usr/lib/x86_64-linux-gnu/'
pythonincpath = '/usr/include/python3.%s'%subversion
boost_libs=['boost_python3%s'%subversion,'boost_numpy3%s'%subversion,'boost_random']
boost_libs='boost_python3%s'%subversion
boost_prefix=['/usr/include','/usr/lib/x86_64-linux-gnu/']
werror=0

mpi = 'OPENMPI'
mpi_prefix = '/usr/lib/openmpi'
mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']
parmetis = True


umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

silo = True
silo_libs = ['siloh5', 'hdf5']
silo_prefix=[ '/usr/include' , '/usr/lib/x86_64-linux-gnu/hdf5/serial', '/usr/lib/x86_64-linux-gnu']

netcdf = 4

netcdf = True
lapack_prefix = ['/usr/include/atlas', '/usr/lib/atlas-base']

from site_init import getdebbuildflags
# Now we add the debian build flags
debstuff = getdebbuildflags()
if len(debstuff) > 0:
  print("Building with the following additional flags from debian: "+str(debstuff))
for i in debstuff:
  k=i[0]
  v=i[1]
  try:
    exec(k+"+=' "+v+"'")
  except NameError:   
    exec(k+"='"+v+"'")


cxx_extra=' -fdiagnostics-color=always'
cxx_extra+=' -Wno-format-truncation'



