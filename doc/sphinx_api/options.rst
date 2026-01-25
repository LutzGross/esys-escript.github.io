=========================
Build Options Reference
=========================

This page provides a complete reference for all build configuration options available in esys-escript.


Option Files
============

The directory ``scons/templates`` contains template configuration files for escript
for various platforms and configurations.


Configuring the Build
---------------------

To build escript on a specific machine create a file named::

    <sourcedir>/scons/<hostname>_options.py

where ``<sourcedir>`` is the escript source directory and ``<hostname>`` is the
machine's short hostname.

If you find a template file whose name matches what you are running, you can import
that file from your new options file then customize to your needs. For example, if you
are running Ubuntu and would like to build with MPI support, you can insert the following
in your options file:

.. code-block:: python

   from templates.ubuntu_options import *
   mpi = 'OPENMPI'

If you can't find a matching template file you can either import one that comes close
or start from scratch and set the options as required.

Delete all variables that are not option variables. In some cases, variables that are
introduced in option files and that are not option variables can trigger an error in scons.


Prefixes
========

There are two ways to specify where to find dependency headers and libraries
(via the ``<dependency>_prefix`` option):

1. If your installation follows the general scheme where headers are located
   in ``<prefix>/include32`` or ``<prefix>/include64``, and libraries in
   ``<prefix>/lib32`` or ``<prefix>/lib64`` then it is sufficient to specify
   this prefix, e.g. ``boost_prefix='/usr'``

2. Otherwise, provide a list with two elements, where the first one is the
   include path, and the second the library path, e.g.

   .. code-block:: python

      boost_prefix=['/usr/include/boost1_48', '/usr/lib']

All ``<dependency>_prefix`` settings default to ``'/usr'``.


Options Reference
=================

The following is an exhaustive list of escript options you can set.
The final setting can be checked by running:

.. code-block:: bash

   scons -h


Core Options
------------

``escript_opts_version``
   The options file version. SCons will refuse to build if there have been changes
   to the set of variables and your file has not been updated. **This setting is mandatory.**

   Current version: ``203``

``prefix``
   Installation prefix - files will be installed in subdirectories underneath this path.

   Default: ``'<sourcedir>'`` (source directory)

``build_dir``
   Top-level directory for intermediate build and test files.

   Default: ``'<sourcedir>/build'``

``verbose``
   Set to ``True`` to print the full compiler/linker command line.

   Default: ``False``


Compiler Options
----------------

``cxx``
   C++ compiler command name or full path.

   Default: auto-detected

``cc_flags``
   Flags to use with the C++ compiler. Do not set this unless you know what you are
   doing - use ``cxx_extra`` to specify additional flags.

   Default: compiler-dependent

``cc_optim``
   Additional compiler (optimization) flags, only applied for non-debug builds.

   Example: ``'-O3 -march=native'``

   Default: compiler-dependent

``cc_debug``
   Additional compiler flags only applied for debug builds.

   Example: ``'-g3 -fno-omit-frame-pointer -D_GLIBCXX_DEBUG'``

   Default: compiler-dependent

``cxx_extra``
   Additional flags to add to the C++ compiler.

   Example: ``'-Wextra -Wno-unused-parameter -I/opt/local/include'``

   Default: ``''`` (empty)

``ld_extra``
   Additional flags to add to the linker.

   Default: ``''`` (empty)

``werror``
   Whether to treat compiler warnings as errors.

   Default: ``True``

``debug``
   Whether to build a debug version (applying ``cc_debug`` flags).

   Default: ``False``


CUDA Options
------------

``nvcc``
   Path to CUDA compiler.

   Default: auto-detected

``nvccflags``
   Flags for CUDA compiler.

   Example: ``'-arch=sm_35 -DBOOST_NOINLINE="__attribute__((noinline))"'``

   Default: ``''`` (empty)

``cuda``
   Whether to add support for GPU-based ripley system matrix (requires nvcc and thrust headers, experimental).

   Default: ``False``

``cuda_prefix``
   Prefix or paths to NVidia CUDA installation.

   Default: ``'/usr/local'``


OpenMP Options
--------------

``openmp``
   Set to ``True`` to add flags that enable OpenMP parallelization.

   Default: ``False``

``omp_flags``
   Additional compiler flags for OpenMP builds.

   Example: ``'-fopenmp'``

   Default: compiler-dependent

``omp_ldflags``
   Additional linker flags for OpenMP builds.

   Example: ``'-fopenmp'``

   Default: compiler-dependent


Boost Options
-------------

``boost_prefix``
   Prefix or paths to boost headers and libraries.

   Example: ``'/usr/local'``

``boost_libs``
   boost-python library/libraries to link against (Python 3).

   Example: ``['boost_python3', 'boost_numpy3']``


MPI Options
-----------

``mpi``
   Flavour of MPI implementation.

   Recognized values:

   * ``'auto'``: Auto-detect MPI implementation from mpi4py (recommended when ``mpi4py=True``)

     - If ``mpi4py=True``: Detects MPI flavour from installed mpi4py package
     - If ``mpi4py=False``: Sets ``mpi='none'``

   * ``'none'``: Disable MPI
   * ``'MPT'``, ``'MPICH'``, ``'MPICH2'``, ``'OPENMPI'``, ``'INTELMPI'``: Explicitly specify MPI implementation

   Default: ``'none'`` (disable MPI)

   .. note::

      When ``mpi4py=True``, the MPI flavour must match the MPI implementation that mpi4py was compiled against.
      Using ``mpi='auto'`` handles this automatically.

``mpi_prefix``
   Prefix or paths to MPI headers and libraries.

   Example: ``'/usr/lib/openmpi'``

``mpi_libs``
   MPI libraries to link against.

   Example: ``['mpi_cxx', 'mpi', 'open-rte', 'open-pal']``

``mpi_no_host``
   This is a workaround for a problem with the open-mpi launcher. In some versions,
   ``"--host localhost"`` would be added to the launch cmd and this would prevent
   multiple processes from being spawned. You should not need to set this option to
   ``True`` unless you are having said problem.

   Default: ``False``

``mpi4py``
   Whether to build escript with the *mpi4py* library. When enabled, allows passing
   custom MPI communicators to domain factory functions via the ``comm`` parameter.

   Default: ``False``


Library Options
---------------

CppUnit
^^^^^^^

``cppunit_prefix``
   Prefix or paths to CppUnit headers and libraries. Only required if you like to run the C++ unit tests.

``cppunit_libs``
   CppUnit library/libraries to link against.

   Default: ``['cppunit']``

HDF5
^^^^

``hdf5``
   Whether to use the HDF5 (serial) library for dump file support.

   Default: ``False``

``hdf5_prefix``
   Prefix or paths to HDF5 headers and libraries.

   Example: ``['/usr/include/hdf5', '/usr/lib']``

``hdf5_libs``
   HDF5 library/libraries to link against.

   Default: ``['hdf5', 'hdf5']``

ParMETIS
^^^^^^^^

``parmetis``
   Whether to use the parMETIS library (only relevant if building finley with MPI).

   Default: ``False``

``parmetis_prefix``
   Prefix or paths to parMETIS headers and libraries.

``parmetis_libs``
   parMETIS library/libraries to link against.

   Default: ``['parmetis', 'metis']``

Intel MKL
^^^^^^^^^

``mkl``
   Whether to add support for the Intel MKL (Math Kernel Library) direct solver.

   Default: ``False``

``mkl_prefix``
   Prefix or paths to MKL headers and libraries.

   Example: ``['/opt/intel/mkl/include', '/opt/intel/mkl/lib/intel64']``

``mkl_libs``
   MKL library/libraries to link against.

   Example: ``['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'pthread']``

UMFPACK
^^^^^^^

``umfpack``
   Whether to add support for the UMFPACK direct solver (requires AMD and BLAS).

   Default: ``False``

``umfpack_prefix``
   Prefix or paths to UMFPACK headers and libraries.

   Example: ``['/usr/include/suitesparse', '/usr/lib']``

``umfpack_libs``
   UMFPACK library/libraries to link against.

   Default: ``['umfpack', 'blas', 'amd']``

MUMPS
^^^^^

``mumps_seq``
   Whether to add support for the sequential MUMPS direct solver.

   Default: ``False``

``mumps_seq_prefix``
   Prefix or paths to sequential MUMPS headers and libraries.

``mumps_seq_libs``
   Sequential MUMPS library/libraries to link against.

   Default: ``['mumps_common', 'dmumps', 'zmumps']``

LAPACK
^^^^^^

``lapack``
   Whether to use BLAS/LAPACK. Note, LAPACK is incompatible with long indices.
   ``'auto'`` tries to detect MKL LAPACK or LAPACKE (modern C interface).

   Default: ``'auto'``

``lapack_prefix``
   Prefix or paths to LAPACK headers and libraries.

   Example: ``['/usr/include', '/usr/lib/x86_64-linux-gnu']``

``lapack_libs``
   LAPACK library/libraries to link against.

   Default: ``'lapacke'`` (modern C interface)

   On Debian/Ubuntu, install with: ``sudo apt-get install liblapacke-dev``

NetCDF
^^^^^^

``netcdf``
   Whether to use NetCDF library for reading/writing NetCDF grid files.
   Value of ``4`` indicates NetCDF version 4.

   Default: ``False``

``netcdf_prefix``
   Prefix or paths to NetCDF headers and libraries.

``netcdf_libs``
   NetCDF library/libraries to link against.

   Default: ``['netcdf_c++4', 'netcdf']``

   On Debian/Ubuntu, install with: ``sudo apt-get install libnetcdf-dev libnetcdf-c++4-dev``

SILO
^^^^

``silo``
   Whether to use LLNL's SILO library for Silo output file support in weipa.

   Default: ``False``

``silo_prefix``
   Prefix or paths to SILO headers and libraries.

``silo_libs``
   SILO library/libraries to link against.

   Default: ``['siloh5', 'hdf5']``


Trilinos Options
----------------

``trilinos``
   Whether to enable support for the Trilinos solver stack.

   Default: ``False``

``build_trilinos``
   Whether to build the trilinos library.

   Recognized values: ``'make'``, ``'always'``, ``'check'``, ``'never'``

   Default: ``'make'``

``trilinos_build``
   The build directory for Trilinos.

``trilinos_prefix``
   Prefix or paths to Trilinos headers and libraries.

``trilinos_libs``
   Trilinos libraries to link against.

   Default: auto-detected


VisIt Options
-------------

``visit``
   Whether to use LLNL's VisIt simulation interface (only version 2 supported).

   Default: ``False``

``visit_prefix``
   Prefix or paths to VisIt's sim2 headers and libraries.

   Example: ``'/opt/visit/2.1.0/linux-intel/libsim/V2'``

``visit_libs``
   Sim2 library/libraries to link against.

   Default: ``['simV2']``


Build Component Options
-----------------------

``domains``
   List of domain families to build.

   Default: ``['finley', 'ripley', 'speckley', 'oxley']``

``paso``
   Whether to build the Paso solver library. Setting this to ``False`` only makes
   sense if you have Trilinos enabled.

   Default: ``True``

``weipa``
   Whether to build the weipa data export library.

   Default: ``True``

``p4est``
   Whether to build escript with the p4est library. This library is required by oxley.

   Default: ``False``

zlib
^^^^

``zlib``
   Whether to enable zlib compression library support. Required by p4est for the oxley domain.

   Default: ``False``

   On Debian/Ubuntu, install with: ``sudo apt-get install zlib1g-dev``

``zlib_prefix``
   Prefix or paths to zlib headers and libraries.

   Default: ``'/usr'``

``zlib_libs``
   zlib library/libraries to link against.

   Default: ``['zlib']``

``sympy``
   Whether to build escript with sympy symbolic module support, if it is available.

   Default: ``False``


Advanced Options
================

.. warning::

   Setting the following options may break your build.

``prelaunch``
   Pre-launch command string.

   Example: ``"EE=$(echo %e|sed -e 's/,/ -x /g')"``

``launcher``
   Launcher command string.

   Example: ``"mpirun -x ${EE} --bynode --bind-to-none --host %h -np %N %b"``

``postlaunch``
   Post-launch command string.

   Default: ``""``

With these three options you can define how to launch programs in your environment.
This is relevant for MPI builds and/or where a batch system or job scheduler is in use.

The content of these options is literally copied into the escript launcher after
applying the following substitutions:

* ``%b`` = executable
* ``%n`` = number of nodes
* ``%p`` = number of processes
* ``%N`` = total number of processes
* ``%t`` = number of threads
* ``%f`` = name of hostfile
* ``%h`` = comma-separated list of hosts
* ``%e`` = comma-separated list of environment variables to export


Other Advanced Options
----------------------

``iknowwhatimdoing``
   Enables code that is non-standard and not recommended for general use.

   Default: ``False``

``insane``
   For testing use only: Skips the sanity check after compilation.

   Default: ``False``

``tools_names``
   Compiler toolset to use.

   Example: ``['intelc']``

   Default: auto-detected

``env_export``
   Additional environmental variables to export to the tools.

   Default: ``[]``

``forcelazy``
   For testing use only, sets the default value for autolazy.

   Default: ``'leave_alone'``

``forcecollres``
   For testing use only, sets the default value for force resolving collective operations.

   Default: ``'leave_alone'``

``sys_libs``
   Extra libraries to link with.

   Default: ``[]``


Python Options
--------------

``pythoncmd``
   Python executable to use for compiling. Must be compatible with the boost python library.

   Default: auto-detected (interpreter executing scons)

``pythonlibname``
   Name of the Python library.

   Example: ``'python3.5m'``

   Default: auto-detected

``pythonlibpath``
   Path to Python library.

   Default: auto-detected

``pythonincpath``
   Path to Python include files.

   Default: auto-detected


Miscellaneous Options
---------------------

``longindices``
   Whether to map ``index_t`` to long (for very large local matrices).

   Default: ``False``

``compressed_files``
   Enable reading compressed binary grids in ripley? (requires boost iostreams)

   Default: ``True``

``compression_libs``
   Compression libraries to link with.

   Default: ``'boost_iostreams'``

``disable_boost_numpy``
   Do not build escript with the boost numpy libraries, even if they are available.

   Default: ``False``

``osx_dependency_fix``
   Whether to apply a dependency fix to libraries (only relevant on OS X).

   Default: ``False``
