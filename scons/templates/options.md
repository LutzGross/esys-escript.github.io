# Option files

The directory [scons/templates](.) contains
template configuration files for escript for various platforms and configurations.

## Configuring build

To build escript on a specific machine create a file named

    <sourcedir>/scons/<hostname>_options.py

where `<sourcedir>` is the escript source directory and `<hostname>` is the
machine's short hostname.
If you find a template file whose name matches what you are
running, you can import that file from your new options file
then customize to your needs. For example, if you are running Ubuntu and would like to build with MPI
support, you can insert the following in your options file:

    from templates.ubuntu_options import *
    mpi = 'OPENMPI'

If you can't find a matching template file you can either 
import one that comes close or start from scratch and set 
the options as required. All recognised options are
explained below.

Delete all variables that are not option variables. In some cases, variables
that are introduced in option files and that are not option variables can trigger
an error in scons.

## Prefixes

There are two ways to specify where to find dependency headers and libraries
(via the `<dependency>_prefix` option):

- If your installation follows the general scheme where headers are located
   in `<prefix>/include32` or `<prefix>/include64`, and 
   libraries in `<prefix>/lib32` or `<prefix>/lib64` then
   it is sufficient to specify this prefix, e.g. `boost_prefix='/usr'`
- Otherwise, provide a list with two elements, where the first one is the
   include path, and the second the library path, e.g.

   boost_prefix=['/usr/include/boost1_48', '/usr/lib']

All `<dependency>_prefix` settings default to '/usr'

## Options

The following is an exhaustive list of escript options you can set.
The final setting can be checked by running:

     scons -h

Each option is followed by a brief explanation.

- `escript_opts_version = 203`:
  The options file version. SCons will refuse to build if there have been
  changes to the set of variables and your file has not been updated.
  This setting is mandatory.

- `prefix = '/usr/local'`:
  Installation prefix - files will be installed in subdirectories underneath
  this path. DEFAULT: '<sourcedir>' (source directory)

- `build_dir = '/tmp/escriptbuild'`:
  Top-level directory for intermediate build and test files.
  DEFAULT: '<sourcedir>/build'

- `verbose = True`:
  Set to True to print the full compiler/linker command line.  DEFAULT: False

- `cxx = 'g++'`:
  C++ compiler command name or full path. DEFAULT: auto-detected

- `cc_flags = ''`:
  Flags to use with the C++ compiler. Do not set this unless you know
  what you are doing - use cxx_extra to specify additional flags.
  DEFAULT: compiler-dependent

- `cc_optim = '-O3 -march=native'`:
  Additional compiler (optimization) flags, only applied for non-debug builds
  DEFAULT: compiler-dependent

- `cc_debug = '-g3 -fno-omit-frame-pointer -D_GLIBCXX_DEBUG'`:
  Additional compiler flags only applied for debug builds
  DEFAULT: compiler-dependent

- `cxx_extra = '-Wextra -Wno-unused-parameter -I/opt/local/include'`:
  Additional flags to add to the C++ compiler. DEFAULT: '' (empty)

- `ld_extra = ''`:
  Additional flags to add to the linker. DEFAULT: '' (empty)

- `nvcc = '/usr/local/bin/nvcc'`:
  Path to CUDA compiler. DEFAULT: auto-detected

- `nvccflags = '-arch=sm_35 -DBOOST_NOINLINE="__attribute__((noinline))"'`:
  Flags for CUDA compiler [new in 202].  DEFAULT: '' (empty)

- `werror = False`:
  Whether to treat compiler warnings as errors. DEFAULT: True

- `debug = False`:
  Whether to build a debug version (applying cc_debug flags)

- `openmp = True`:
  Set to True to add flags that enable OpenMP parallelization
  DEFAULT: False

- `omp_flags = '-fopenmp'`:
  Additional compiler flags for OpenMP builds. DEFAULT: compiler-dependent

- `omp_ldflags = '-fopenmp'`:
  Additional linker flags for OpenMP builds. DEFAULT: compiler-dependent

- `boost_prefix = '/usr/local'`:
  Prefix or paths to boost headers and libraries. See note above.

- `boost_libs = ['boost_python3', 'boost_numpy3']`:
  boost-python library/libraries to link against (Python 3)

- `build_trilinos = 'never'`:
  Whether to build the trilinos library
  Recognized values: 'make', 'always', 'check', 'never'
  DEFAULT: 'make'

- `cppunit_prefix = '/usr/local'`:
  Prefix or paths to CppUnit headers and libraries. See note above.
  Only required if you like to run the C++ unit tests.

- `cppunit_libs = ['cppunit']`:
  CppUnit library/libraries to link against. Only required if you like to run
  the C++ unit tests.

- `mpi = 'OPENMPI'`:
  Flavour of MPI implementation
  Recognized values: 'auto', 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI'
  - 'auto': Auto-detect MPI implementation from mpi4py (recommended when mpi4py=True)
    * If mpi4py=True: Detects MPI flavour from installed mpi4py package
    * If mpi4py=False: Sets mpi='none'
  - 'none': Disable MPI
  - Other values: Explicitly specify MPI implementation
  DEFAULT: 'none' (disable MPI)
  **Note:** When mpi4py=True, the MPI flavour must match the MPI implementation that mpi4py was compiled against. Using mpi='auto' handles this automatically.

- `mpi_prefix = '/usr/lib/openmpi'`:
  Prefix or paths to MPI headers and libraries. See note above about prefixes.

- `mpi_libs = ['mpi_cxx', 'mpi', 'open-rte', 'open-pal']`:
  MPI libraries to link against.

- `mpi_no_host = False`:
  This is a workaround for a problem with the open-mpi launcher. In some 
  versions, `"--host localhost"` would be added to the launch cmd and this
  would prevent multiple processes from being spawned.
  You should not need to set this option to True unless you are having
  said problem.

- `mpi4py = False`:
  Whether to build escript with the *mpi4py* library. When enabled, allows passing custom MPI communicators to domain factory functions via the `comm` parameter.
  DEFAULT: False

- `cuda = True`:
  Whether to add support for GPU-based ripley system matrix (requires nvcc
  and thrust headers, experimental).  DEFAULT: False

- `cuda_prefix = '/usr/local'`:
  Prefix or paths to NVidia CUDA installation. See note above. [new in 202]

- `hdf5 = False`:
  Whether to use the HDF5 (serial) library for dump file support.

- `hdf5_prefix = ['/usr/include/hdf5', '/usr/lib']`:
  Prefix or paths to HDF5 headers and libraries. See note above.

- `hdf5_libs = ['hdf5', 'hdf5']`:
  HDF5 library/libraries to link against

- `parmetis = True`:
  Whether to use the parMETIS library (only relevant if building
  finley with MPI). DEFAULT: False

- `parmetis_prefix = '/usr/local'`:
  Prefix or paths to parMETIS headers and libraries. See note above.

- `parmetis_libs = ['parmetis', 'metis']`:
  parMETIS library/libraries to link against

- `mkl = False`:
  Whether to add support for the Intel MKL (Math Kernel Library) direct solver
  
- `mkl_prefix = ['/opt/intel/mkl/include', '/opt/intel/mkl/lib/intel64']`:
  Prefix or paths to MKL headers and libraries. See note above.

- `mkl_libs = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']`:
  MKL library/libraries to link against

- `umfpack = False`:
  Whether to add support for the UMFPACK direct solver (requires AMD and BLAS)
  DEFAULT: False

- `umfpack_prefix = ['/usr/include/suitesparse', '/usr/lib']`:
  Prefix or paths to UMFPACK headers and libraries. See note above.

- `umfpack_libs = ['umfpack', 'blas', 'amd']`:
  UMFPACK library/libraries to link against

- `mumps_seq = False`:
  Whether to add support for the sequential MUMPS direct solver

- `mumps_seq_prefix = ['/usr/lib']`:
  Prefix or paths to sequential MUMPS headers and libraries. See note above.

- `mumps_seq_libs = ['mumps_common', 'dmumps', 'zmumps']`:
  Sequential MUMPS library/libraries to link against


- `p4est = False`:
  Whether to build escript with the p4est library. This library is required by oxley.

- `zlib = False`:
  Whether to enable zlib compression library support. Required by p4est for the oxley domain.
  On Debian/Ubuntu, install with: `sudo apt-get install zlib1g-dev`

- `zlib_prefix = '/usr'`:
  Prefix or paths to zlib headers and libraries.

- `zlib_libs = ['zlib']`:
  zlib library/libraries to link against.

- `sympy = False`:
  Whether to build escript with sympy symbolic module support, if it is available
  DEFAULT: False

- `lapack = auto`:
  Whether to use BLAS/LAPACK. Note, LAPACK is incompatible with long indices.
  `auto` tries to detect MKL LAPACK or LAPACKE (modern C interface).
  DEFAULT: 'auto'

- `lapack_prefix = ['/usr/include', '/usr/lib/x86_64-linux-gnu']`:
  Prefix or paths to LAPACK headers and libraries. See note above.

- `lapack_libs = ['lapacke']`:
  LAPACK library/libraries to link against. Default is 'lapacke' (modern C interface).
  On Debian/Ubuntu, install with: `sudo apt-get install liblapacke-dev`

- `silo = False`:
  Whether to use LLNL's SILO library for Silo output file support in weipa

- `silo_prefix = '/usr/local'`:
  Prefix or paths to SILO headers and libraries. See note above.

- `silo_libs = ['siloh5', 'hdf5']`:
  SILO library/libraries to link against

- `netcdf = 4`:
  Whether to use NetCDF library for reading/writing NetCDF grid files.
  Value of 4 indicates NetCDF version 4.
  DEFAULT: False

- `netcdf_prefix = '/usr/local'`:
  Prefix or paths to NetCDF headers and libraries. See note above.

- `netcdf_libs = ['netcdf_c++4', 'netcdf']`:
  NetCDF library/libraries to link against.
  On Debian/Ubuntu, install with: `sudo apt-get install libnetcdf-dev libnetcdf-c++4-dev`

- `trilinos = True`:
  Whether to enable support for the Trilinos solver stack. [new in 203]
  DEFAULT: False

- `trilinos_build = '/usr/trilinos_build'`:
  The build directory for Trilinos

- `trilinos_prefix = '/usr/local'`:
  Prefix or paths to Trilinos headers and libraries. See note above.

- `trilinos_libs = []`:
  Trilinos libraries to link against. DEFAULT: auto-detected

- `visit = False`:
  Whether to use LLNL's VisIt simulation interface (only version 2 supported).

- `visit_prefix = '/opt/visit/2.1.0/linux-intel/libsim/V2'`:
  Prefix or paths to VisIt's sim2 headers and libraries. See note above.

- `visit_libs = ['simV2']`:
  Sim2 library/libraries to link against

- `domains = ['finley', 'ripley', 'speckley', 'oxley']`:
  List of domain families to build.

- `paso = True`:
  Whether to build the Paso solver library. Setting this to False only makes
  sense if you have Trilinos enabled. DEFAULT: True

- `weipa = True`:
  Whether to build the weipa data export library. DEFAULT: True

### Advanced Options

Setting the following options may break your build.

- `prelaunch = "EE=$(echo %e|sed -e 's/,/ -x /g')"`

- `launcher = "mpirun -x ${EE} --bynode --bind-to-none --host %h -np %N %b"`
- `postlaunch = ""`
  With these three options you can define how to launch programs in your
  environment. This is relevant for MPI builds and/or where a batch system
  or job scheduler is in use. 
  The content of these options is literally copied into the escript launcher
  after applying the following substitutions:
  
     - %b = executable, %n = number of nodes, %p = number of processes,
     - %N = total number of processes, %t = number of threads,
     - %f = name of hostfile, %h = comma-separated list of hosts,
     - %e = comma-separated list of environment variables to export

- `iknowwhatimdoing = True`:
  Enables code that is non-standard and not recommended for general use.

- `insane = True`:
  For testing use only: Skips the sanity check after compilation
  DEFAULT: 'False'

- `tools_names = ['intelc']`:
  compiler toolset to use. DEFAULT: auto-detected

- `env_export = []`:
  Additional environmental variables to export to the tools

- `forcelazy = 'on'`:
  For testing use only, sets the default value for autolazy.
  DEFAULT: 'leave_alone'

- `forcecollres = 'on'`:
  For testing use only, sets the default value for force resolving collective
  operations.  DEFAULT: 'leave_alone'

- `sys_libs = []`:
  Extra libraries to link with

- `pythoncmd = '/usr/bin/python3'`:
  Python executable to use for compiling. Must be compatible with the
  boost python library
  DEFAULT: auto-detected (interpreter executing scons)

- `pythonlibname = 'python3.5m'`:
  Name of the Python library. DEFAULT: auto-detected.

- `pythonlibpath = '/usr/lib'`:
  Path to Python library. DEFAULT: auto-detected.

- `pythonincpath = '/usr/include/python3.5'`:
  Path to Python include files. DEFAULT: auto-detected.

- `longindices = false`:
  Whether to map index_t to long (for very large local matrices) [new in 202]
  DEFAULT: False

- `compressed_files = False`:
  Enable reading compressed binary grids in ripley? (requires boost iostreams)
  DEFAULT: True

- `compression_libs = ['boost_iostreams-mt']`:
  Compression libraries to link with. DEFAULT: 'boost_iostreams'

- `disable_boost_numpy = False`:
  Do not build escript with the boost numpy libraries, even if they are available.
  DEFAULT: False

- `osx_dependency_fix = True`:
  Whether to apply a dependency fix to libraries (only relevant on OS X).
  DEFAULT: False

