
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


# locations of include files for python
### python_path = '/data/raid2/toolspp4/python/2.4.3/gcc-3.3.6/include/python2.4'
### python_lib_path = '/data/raid2/toolspp4/python/2.4.3/gcc-3.3.6/lib'
### python_lib = 'python2.4'

# locations of libraries for boost
### boost_path = '/data/raid2/toolspp4/boost/1.33.1/python-2.4.3/gcc-3.3.6/include'
### boost_lib_path = '/data/raid2/toolspp4/boost/1.33.1/python-2.4.3/gcc-3.3.6/lib'
### boost_lib = 'boost_python-mt'

# locations of doc building executables
### doxygen_path = '/data/raid2/toolspp4/doxygen/1.4.6/gcc-3.3.6/bin'
### epydoc_path = '/raid2/tools/epydoc/2.1/python-2.3.4/bin'

# locations of netcdf
useNetCDF = 'yes'
netCDF_path = "/usr/include/netcdf-3"
netCDF_lib_path = "/usr/lib"
netCDF_libs = [ 'netcdf_c++', 'netcdf']

### mpi_path = '/usr/include'
### mpi_lib_path = '/usr/lib'
### mpi_libs = [ 'mpi' ]
### mpi_flavour = 'MPICH'

### omp_flags = '-openmp -openmp_report2 '
### omp_flags_debug = '-openmp -openmp_report0'

# c flags to use
### cc_flags  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -c99 -w1 -wd161 -fpic -ivdep-parallel"
### cc_flags_debug  = '-g -O0 -c99 -w1 -wd161 -fpic'

# c++ flags to use
### cxx_flags = '-ansi -wd161 -DMPI_NO_CPPBIND'
### cxx_flags_debug = '-ansi -wd161 -DDOASSERT -DDOPROF -DMPI_NO_CPPBIND'

# c and c++ flags for MPI compilation
# c flags to use
### cc_flags_MPI  = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -c99 -w1 -fpic -wd161 -DPASO_MPI -ivdep-parallel"
### cc_flags_debug_MPI  = '-g -O0 -c99 -w1 -fpic -wd161 -DPASO_MPI'

# c++ flags to use
### cxx_flags_MPI = '-ansi -wd1563 -wd161 -DMPI_NO_CPPBIND'
### cxx_flags_debug_MPI = '-ansi -DDOASSERT -DDOPROF -wd1563 -wd161 -DMPI_NO_CPPBIND'

# system specific libraries to link with
### sys_libs = ['guide', 'irc']
