##############################################################################
#
# GitHub Actions Ubuntu build configuration for esys-escript
# Minimal build for documentation generation (no Trilinos)
#
##############################################################################

# MANDATORY: Options file version
escript_opts_version = 203

# Installation and build directories
build_dir = 'build'

# Compiler settings
openmp = True
werror = False

# Boost configuration
boost_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu/']
# boost_libs will be set dynamically based on Python version

# Direct solvers
umfpack = True
umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']
umfpack_libs = ['umfpack', 'blas', 'amd']

# File I/O libraries
hdf5 = True
hdf5_prefix = ['/usr/include/hdf5/serial', '/usr/lib/x86_64-linux-gnu/hdf5/serial']
hdf5_libs = ['hdf5_serial', 'hdf5_serial_cpp']

silo = True
silo_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']
silo_libs = ['siloh5', 'hdf5_serial', 'hdf5_serial_cpp']

# Compression support
compressed_files = True
# compression_libs will be set based on boost version

# Python configuration
pythoncmd = 'python3'

# Solver configuration
paso = True
trilinos = False  # Disabled for faster build
build_trilinos = 'never'  # Don't build Trilinos from source

# Domain families to build
domains = ['finley', 'ripley', 'speckley' ] # + [ 'oxley']

# Export library
weipa = True

# Disable MPI for GitHub Actions
mpi = 'none'

# No sympy or other optional features
sympy = False

# Disable netCDF - not needed for docs and causes header issues if not properly configured
netcdf = False
