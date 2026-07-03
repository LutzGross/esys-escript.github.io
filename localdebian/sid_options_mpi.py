# This is a template configuration file for escript on Debian/GNU Linux.
# Refer to README_FIRST for usage instructions.

from debian.sid_options import *

# mpi = 'OPENMPI'
mpi_include = '/usr/lib/' + arch + '/' + mpi.lower() + '/include'
mpi_lib = '/usr/lib/' + arch + '/' + mpi.lower() + '/lib'
mpi_prefix = [mpi_include,mpi_lib]

cxx_extra=' -Wno-stringop-truncation'

