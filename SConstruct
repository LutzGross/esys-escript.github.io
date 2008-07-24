#         Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php

# TODO How to modify CCFLAGS for only one file? Look for example with Program(...)

EnsureSConsVersion(0,96,91)
EnsurePythonVersion(2,3)

import sys, os, re, socket

# Add our extensions
if os.path.isdir('scons'): sys.path.append('scons')
import scons_extensions

# The string python2.4 or python2.5
python_version = 'python%s.%s' % (sys.version_info[0], sys.version_info[1])

# MS Windows support, many thanks to PH
IS_WINDOWS_PLATFORM = (os.name== "nt")

prefix = ARGUMENTS.get('prefix', Dir('#.').abspath)

# Read configuration options from file scons/<hostname>_options.py
hostname = re.sub("[^0-9a-zA-Z]", "_", socket.gethostname().split('.')[0])
tmp = os.path.join("scons",hostname+"_options.py")
options_file = ARGUMENTS.get('options_file', tmp)
if not os.path.isfile(options_file): options_file = False

# Load options file and command-line arguments
opts = Options(options_file, ARGUMENTS)

############ Load build options ################################

opts.AddOptions(
# Where to install esys stuff
  ('prefix', 'where everything will be installed',                       Dir('#.').abspath),
  ('incinstall', 'where the esys headers will be installed',             Dir('#.').abspath+'/include'),
  ('libinstall', 'where the esys libraries will be installed',           os.path.join(prefix,"lib")),
  ('pyinstall', 'where the esys python modules will be installed',       os.path.join(prefix,"esys")),
# Compilation options
  BoolOption('dodebug', 'For backwards compatibility', 'no'),
  BoolOption('usedebug', 'Do you want a debug build?', 'no'),
  BoolOption('usevtk', 'Do you want to use VTK?', 'yes'),
  ('options_file', "File of paths/options. Default: scons/<hostname>_options.py", options_file),
  # The strings -DDEFAULT_ get replaced by scons/<hostname>_options.py or by defaults below
  ('cc_flags', 'C compiler flags to use', '-DEFAULT_1'),
  ('cc_optim', 'C compiler optimization flags to use', '-DEFAULT_2'),
  ('cc_debug', 'C compiler debug flags to use', '-DEFAULT_3'),
  ('omp_optim', 'OpenMP compiler flags to use (Release build)', '-DEFAULT_4'),
  ('omp_debug', 'OpenMP compiler flags to use (Debug build)', '-DEFAULT_5'),
  ('sys_libs', 'System libraries to link with', []),
  ('ar_flags', 'Static library archiver flags to use', ''),
  BoolOption('useopenmp', 'Compile parallel version using OpenMP', 'yes'),
# Python
  ('python_path', 'Path to Python includes', '/usr/include/'+python_version),
  ('python_lib_path', 'Path to Python libs', '/usr/lib'),
  ('python_libs', 'Python libraries to link with', [python_version]),
  ('python_cmd', 'Python command', 'python'),
# Boost
  ('boost_path', 'Path to Boost includes', "/usr/include"),
  ('boost_lib_path', 'Path to Boost libs', '/usr/lib'),
  ('boost_libs', 'Boost libraries to link with', ['boost_python']),
# NetCDF
  BoolOption('usenetcdf', 'switch on/off the usage of netCDF', 'yes'),
  ('netCDF_path', 'Path to netCDF includes', '/usr/include'),
  ('netCDF_lib_path', 'Path to netCDF libs', '/usr/lib'),
  ('netCDF_libs', 'netCDF C++ libraries to link with', ['netcdf_c++', 'netcdf']),
# MPI
  BoolOption('useMPI', 'For backwards compatibility', 'no'),
  BoolOption('usempi', 'Compile parallel version using MPI', 'no'),
  ('MPICH_IGNORE_CXX_SEEK', 'name of macro to ignore MPI settings of C++ SEEK macro (for MPICH)' , 'MPICH_IGNORE_CXX_SEEK'),
  ('mpi_path', 'Path to MPI includes', '/usr/include'),
  ('mpi_run', 'mpirun name' , 'mpiexec -np 1'),
  ('mpi_lib_path', 'Path to MPI libs (needs to be added to the LD_LIBRARY_PATH)', '/usr/lib'),
  ('mpi_libs', 'MPI libraries to link with (needs to be shared!)', ['mpich' , 'pthread', 'rt']),
# ParMETIS
  BoolOption('useparmetis', 'Compile parallel version using ParMETIS', 'yes'),
  ('parmetis_path', 'Path to ParMETIS includes', '/usr/include'),
  ('parmetis_lib_path', 'Path to ParMETIS library', '/usr/lib'),
  ('parmetis_libs', 'ParMETIS library to link with', []),
# PAPI
  BoolOption('usepapi', 'switch on/off the usage of PAPI', 'no'),
  ('papi_path', 'Path to PAPI includes', '/usr/include'),
  ('papi_lib_path', 'Path to PAPI libs', '/usr/lib'),
  ('papi_libs', 'PAPI libraries to link with', ['papi']),
  BoolOption('papi_instrument_solver', 'use PAPI in Solver.c to instrument each iteration of the solver', False),
# MKL
  BoolOption('usemkl', 'switch on/off the usage of MKL', 'yes'),
  ('mkl_path', 'Path to MKL includes', '/opt/intel/mkl80.019/include'),
  ('mkl_lib_path', 'Path to MKL libs', '/opt/intel/mkl80.019/lib/64'),
  ('mkl_libs', 'MKL libraries to link with', ['mkl_solver', 'mkl_lapack', 'mkl_em64t']),
# UMFPACK
  BoolOption('useumfpack', 'switch on/off the usage of UMFPACK', 'yes'),
  ('ufc_path', 'Path to UFconfig includes', '/usr/include/suitesparse'),
  ('umf_path', 'Path to UMFPACK includes', '/usr/include/suitesparse'),
  ('umf_lib_path', 'Path to UMFPACK libs', '/usr/lib'),
  ('umf_libs', 'UMFPACK libraries to link with', ['umfpack']),
# AMD (used by UMFPACK)
  ('amd_path', 'Path to AMD includes', '/usr/include/suitesparse'),
  ('amd_lib_path', 'Path to AMD libs', '/usr/lib'),
  ('amd_libs', 'AMD libraries to link with', ['amd']),
# BLAS (used by UMFPACK)
  ('blas_path', 'Path to BLAS includes', '/usr/include/suitesparse'),
  ('blas_lib_path', 'Path to BLAS libs', '/usr/lib'),
  ('blas_libs', 'BLAS libraries to link with', ['blas'])
)

############ Specify which compilers to use ####################

# intelc uses regular expressions improperly and emits a warning about
# failing to find the compilers.  This warning can be safely ignored.

if IS_WINDOWS_PLATFORM:
      env = Environment(tools = ['default', 'msvc'], options = opts)
else:
   if socket.gethostname().split('.')[0] == 'service0':
      env = Environment(tools = ['default', 'intelc'], options = opts)
   elif os.uname()[4]=='ia64':
      env = Environment(tools = ['default', 'intelc'], options = opts)
      if env['CXX'] == 'icpc':
         env['LINK'] = env['CXX'] # version >=9 of intel c++ compiler requires use of icpc to link in C++ runtimes (icc does not)
   else:
      env = Environment(tools = ['default'], options = opts)
Help(opts.GenerateHelpText(env))

############ Fill in compiler options if not set above #########

# Backwards compatibility: allow dodebug=yes and useMPI=yes
if env['dodebug']: env['usedebug'] = 1
if env['useMPI']: env['usempi'] = 1

# Default compiler options (override allowed in hostname_options.py, but should not be necessary)
# For both C and C++ you get: cc_flags and either the optim flags or debug flags

if env["CC"] == "icc":
  # Intel compilers
  cc_flags		= "-fPIC -ansi -wd161 -w1 -vec-report0 -DBLOCKTIMER -DCORE_ID1"
  cc_optim		= "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias"
  cc_debug		= "-g -O0 -UDOASSERT -DDOPROF -DBOUNDS_CHECK"
  omp_optim		= "-openmp -openmp_report0"
  omp_debug		= "-openmp -openmp_report0"
elif env["CC"] == "gcc":
  # GNU C on any system
  cc_flags		= "-fPIC -ansi -ffast-math -Wno-unknown-pragmas -pedantic-errors -Wno-long-long -DBLOCKTIMER"
  cc_optim		= "-O3"
  cc_debug		= "-g -O0 -UDOASSERT -DDOPROF -DBOUNDS_CHECK"
  omp_optim		= ""
  omp_debug		= ""
elif env["CC"] == "cl":
  # Microsoft Visual C on Windows
  cc_flags		= "/FD /EHsc /GR /wd4068 -D_USE_MATH_DEFINES -DDLL_NETCDF"
  cc_optim		= "/O2 /Op /MT /W3"
  cc_debug		= "/Od /RTC1 /MTd /ZI -DBOUNDS_CHECK"
  omp_optim		= ""
  omp_debug		= ""

# If not specified in hostname_options.py then set them here
if env["cc_flags"]	== "-DEFAULT_1": env['cc_flags'] = cc_flags
if env["cc_optim"]	== "-DEFAULT_2": env['cc_optim'] = cc_optim
if env["cc_debug"]	== "-DEFAULT_3": env['cc_debug'] = cc_debug
if env["omp_optim"]	== "-DEFAULT_4": env['omp_optim'] = omp_optim
if env["omp_debug"]	== "-DEFAULT_5": env['omp_debug'] = omp_debug

# OpenMP is disabled if useopenmp=no or both variables omp_optim and omp_debug are empty
if not env["useopenmp"]:
  env['omp_optim'] = ""
  env['omp_debug'] = ""

if env['omp_optim'] == "" and env['omp_debug'] == "": env["useopenmp"] = 0

############ Copy environment variables into scons env #########

try: env['ENV']['OMP_NUM_THREADS'] = os.environ['OMP_NUM_THREADS']
except KeyError: env['ENV']['OMP_NUM_THREADS'] = 1

try: env['ENV']['PATH'] = os.environ['PATH']
except KeyError: pass

try: env['ENV']['PYTHONPATH'] = os.environ['PYTHONPATH']
except KeyError: pass

try: env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']
except KeyError: pass

try: env['ENV']['LIBRARY_PATH'] = os.environ['LIBRARY_PATH']
except KeyError: pass

try: env['ENV']['DISPLAY'] = os.environ['DISPLAY']
except KeyError: pass

try: env['ENV']['XAUTHORITY'] = os.environ['XAUTHORITY']
except KeyError: pass

try: env['ENV']['HOME'] = os.environ['HOME']
except KeyError: pass

# Configure for test suite
env.PrependENVPath('PYTHONPATH', prefix)
env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

############ Set up paths for Configure() ######################

# Make a copy of an environment
# Use env.Clone if available, but fall back on env.Copy for older version of scons
def clone_env(env):
  if 'Clone' in dir(env): return env.Clone()	# scons-0.98
  else:                   return env.Copy()	# scons-0.96

# Add cc option -I<Escript>/trunk/include
env.Append(CPPPATH		= [Dir('include')])

# Add cc option -L<Escript>/trunk/lib
env.Append(LIBPATH		= [Dir('lib')])

env.Append(CPPDEFINES = ['ESCRIPT_EXPORTS', 'FINLEY_EXPORTS'])

# Get the global Subversion revision number for getVersion() method
try:
   global_revision = os.popen("svnversion -n .").read()
   global_revision = re.sub(":.*", "", global_revision)
   global_revision = re.sub("[^0-9]", "", global_revision)
except:
   global_revision="-1"
if global_revision == "": global_revision="-2"
env.Append(CPPDEFINES = ["SVN_VERSION="+global_revision])

# Create a Configure() environment only for checking existence of
# libraries and headers.  Later we throw it away then build the real
# environment.
conf = Configure(clone_env(env))

############ numarray (required) ###############################

try: from numarray import identity
except ImportError: sys.exit(1)

############ python libraries (required) #######################

conf.env.Append(CPPPATH		= [env['python_path']])
conf.env.Append(LIBPATH		= [env['python_lib_path']])
conf.env.Append(LIBS		= [env['python_libs']])

if not conf.CheckCHeader('Python.h'): sys.exit(1)
if not conf.CheckFunc('Py_Main'): sys.exit(1)

# Add python libraries to environment env
env.AppendUnique(CPPPATH = [env['python_path']])
env.AppendUnique(LIBPATH = [env['python_lib_path']])
env.Append(LIBS = [env['python_libs']])

############ boost (required) ##################################

conf.env.Append(CPPPATH		= [env['boost_path']])
conf.env.Append(LIBPATH		= [env['boost_lib_path']])
conf.env.Append(LIBS		= [env['boost_libs']])

if not conf.CheckCXXHeader('boost/python.hpp'): sys.exit(1)
if not conf.CheckFunc('PyObject_SetAttr'): sys.exit(1)

# Add boost to environment env
env.AppendUnique(CPPPATH = [env['boost_path']])
env.AppendUnique(LIBPATH = [env['boost_lib_path']])
env.Append(LIBS = [env['boost_libs']])

############ VTK (optional) ####################################

# You must set up your PYTHONPATH before calling scons

if env['usevtk']:
  try:
    import vtk
    env['usevtk'] = 1
  except ImportError:
    env['usevtk'] = 0

# Add VTK to environment env
if env['usevtk']:
  env.Append(CPPDEFINES = ['USE_VTK'])

############ NetCDF (optional) #################################

# Start a new configure environment that reflects what we've already found
conf.Finish()
conf = Configure(clone_env(env))

if env['usenetcdf']:
  conf.env.Append(CPPPATH	= [env['netCDF_path']])
  conf.env.Append(LIBPATH	= [env['netCDF_lib_path']])
  conf.env.Append(LIBS		= [env['netCDF_libs']])

if env['usenetcdf'] and not conf.CheckCHeader('netcdf.h'): env['usenetcdf'] = 0
if env['usenetcdf'] and not conf.CheckFunc('nc_open'): env['usenetcdf'] = 0

# Add NetCDF to environment env
if env['usenetcdf']:
  env.AppendUnique(CPPPATH = [env['netCDF_path']])
  env.AppendUnique(LIBPATH = [env['netCDF_lib_path']])
  env.Append(LIBS = [env['netCDF_libs']])
  env.Append(CPPDEFINES = ['USE_NETCDF'])

############ PAPI (optional) ###################################

# Start a new configure environment that reflects what we've already found
conf.Finish()
conf = Configure(clone_env(env))

if env['usepapi']:
  conf.env.Append(CPPPATH	= [env['papi_path']])
  conf.env.Append(LIBPATH	= [env['papi_lib_path']])
  conf.env.Append(LIBS		= [env['papi_libs']])

if env['usepapi'] and not conf.CheckCHeader('papi.h'): env['usepapi'] = 0
if env['usepapi'] and not conf.CheckFunc('PAPI_start_counters'): env['usepapi'] = 0

# Add PAPI to environment env
if env['usepapi']:
  env.AppendUnique(CPPPATH = [env['papi_path']])
  env.AppendUnique(LIBPATH = [env['papi_lib_path']])
  env.Append(LIBS = [env['papi_libs']])
  env.Append(CPPDEFINES = ['BLOCKPAPI'])

############ MKL (optional) ####################################

# Start a new configure environment that reflects what we've already found
conf.Finish()
conf = Configure(clone_env(env))

if env['usemkl']:
  conf.env.Append(CPPPATH	= [env['mkl_path']])
  conf.env.Append(LIBPATH	= [env['mkl_lib_path']])
  conf.env.Append(LIBS		= [env['mkl_libs']])

if env['usemkl'] and not conf.CheckCHeader('mkl_solver.h'): env['usemkl'] = 0
if env['usemkl'] and not conf.CheckFunc('PARDISO'): env['usemkl'] = 0

# Add MKL to environment env
if env['usemkl']:
  print "ksteube using MKL"
  env.AppendUnique(CPPPATH = [env['mkl_path']])
  env.AppendUnique(LIBPATH = [env['mkl_lib_path']])
  env.Append(LIBS = [env['mkl_libs']])
  env.Append(CPPDEFINES = ['MKL'])

############ UMFPACK (optional) ################################

# Start a new configure environment that reflects what we've already found
conf.Finish()
conf = Configure(clone_env(env))

if env['useumfpack']:
  conf.env.Append(CPPPATH	= [env['ufc_path']])
  conf.env.Append(CPPPATH	= [env['umf_path']])
  conf.env.Append(LIBPATH	= [env['umf_lib_path']])
  conf.env.Append(LIBS		= [env['umf_libs']])
  conf.env.Append(CPPPATH	= [env['amd_path']])
  conf.env.Append(LIBPATH	= [env['amd_lib_path']])
  conf.env.Append(LIBS		= [env['amd_libs']])
  conf.env.Append(CPPPATH	= [env['blas_path']])
  conf.env.Append(LIBPATH	= [env['blas_lib_path']])
  conf.env.Append(LIBS		= [env['blas_libs']])

if env['useumfpack'] and not conf.CheckCHeader('umfpack.h'): env['useumfpack'] = 0
if env['useumfpack'] and not conf.CheckFunc('umfpack_di_symbolic'): env['useumfpack'] = 0

# Add UMFPACK to environment env
if env['useumfpack']:
  env.AppendUnique(CPPPATH = [env['ufc_path']])
  env.AppendUnique(CPPPATH = [env['umf_path']])
  env.AppendUnique(LIBPATH = [env['umf_lib_path']])
  env.Append(LIBS = [env['umf_libs']])
  env.Append(CPPDEFINES = ['UMFPACK'])

  env.AppendUnique(CPPPATH = [env['amd_path']])
  env.AppendUnique(LIBPATH = [env['amd_lib_path']])
  env.Append(LIBS = [env['amd_libs']])

  env.AppendUnique(CPPPATH = [env['blas_path']])
  env.AppendUnique(LIBPATH = [env['blas_lib_path']])
  env.Append(LIBS = [env['blas_libs']])

############ MPI (optional) ####################################

# Start a new configure environment that reflects what we've already found
conf.Finish()
conf = Configure(clone_env(env))

if env['usempi']:
  conf.env.Append(CPPPATH	= [env['mpi_path']])
  conf.env.Append(LIBPATH	= [env['mpi_lib_path']])
  conf.env.Append(LIBS		= [env['mpi_libs']])

if env['usempi'] and not conf.CheckCHeader('mpi.h'): env['usempi'] = 0
if env['usempi'] and not conf.CheckFunc('MPI_Init'): env['usempi'] = 0

# Create a modified environment for MPI programs
env_mpi = clone_env(env)

# Add MPI to environment env_mpi
if env_mpi['usempi']:
  env_mpi.AppendUnique(CPPPATH = [env['mpi_path']])
  env_mpi.AppendUnique(LIBPATH = [env['mpi_lib_path']])
  env_mpi.Append(LIBS = [env['mpi_libs']])
  env_mpi.Append(CPPDEFINES = ['PASO_MPI', 'MPI_NO_CPPBIND', env_mpi['MPICH_IGNORE_CXX_SEEK']])

############ ParMETIS (optional) ###############################

# Start a new configure environment that reflects what we've already found
conf.Finish()
conf = Configure(clone_env(env_mpi))

if not env['usempi']: env['useparmetis'] = 0

if env['useparmetis']:
  conf.env.Append(CPPPATH	= [env['parmetis_path']])
  conf.env.Append(LIBPATH	= [env['parmetis_lib_path']])
  conf.env.Append(LIBS		= [env['parmetis_libs']])

if env['useparmetis'] and not conf.CheckCHeader('parmetis.h'): env['useparmetis'] = 0
if env['useparmetis'] and not conf.CheckFunc('ParMETIS_V3_PartGeomKway'): env['useparmetis'] = 0

# Add ParMETIS to environment env_mpi
if env['useparmetis']:
  env_mpi.AppendUnique(CPPPATH = [env['parmetis_path']])
  env_mpi.AppendUnique(LIBPATH = [env['parmetis_lib_path']])
  env_mpi.Append(LIBS = [env['parmetis_libs']])
  env_mpi.Append(CPPDEFINES = ['USE_PARMETIS'])

############ Configure finished, now finish environment ########

conf.Finish()

# Enable debug by choosing either cc_debug or cc_optim
if env['usedebug']:
  env.Append(CCFLAGS		= env['cc_debug'])
  env.Append(CCFLAGS		= env['omp_debug'])
else:
  env.Append(CCFLAGS		= env['cc_optim'])
  env.Append(CCFLAGS		= env['omp_optim'])

# Always use cc_flags
env.Append(CCFLAGS		= env['cc_flags'])

# MS Windows
if IS_WINDOWS_PLATFORM:
  env.PrependENVPath('PATH',	[env['boost_lib_path']])
  env.PrependENVPath('PATH',	[env['libinstall']])
  if env['usenetcdf']:
    env.PrependENVPath('PATH',	[env['netCDF_lib_path']])

env.Append(ARFLAGS = env['ar_flags'])

############ Summarize our environment #########################

print ""
print "Summary of configuration (see ./config.log for information)"
print "	Using python libraries"
print "	Using numarray"
print "	Using boost"
if env['usenetcdf']: print "	Using NetCDF"
else: print "	Not using NetCDF"
if env['usevtk']: print "	Using VTK"
else: print "	Not using VTK"
if env['usemkl']: print "	Using MKL"
else: print "	Not using MKL"
if env['useumfpack']: print "	Using UMFPACK"
else: print "	Not using UMFPACK"
if env['useopenmp']: print "	Using OpenMP"
else: print "	Not using OpenMP"
if env['usempi']: print "	Using MPI"
else: print "	Not using MPI"
if env['useparmetis']: print "	Using ParMETIS"
else: print "	Not using ParMETIS (requires MPI)"
if env['usepapi']: print "	Using PAPI"
else: print "	Not using PAPI"
if env['usedebug']: print "	Compiling for debug"
else: print "	Not compiling for debug"
print "	Installing in", prefix
print ""

############ Add some custom builders ##########################

py_builder = Builder(action = scons_extensions.build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed', src_suffix=env['PROGSUFFIX'], single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = scons_extensions.runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

############ Build the desired subdirectories ##################

Export(["env", "env_mpi"])

env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
env.SConscript(dirs = ['escript/src'], build_dir='build/$PLATFORM/escript', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['finley/src'], build_dir='build/$PLATFORM/finley', duplicate=0)
env.SConscript(dirs = ['modellib/py_src'], build_dir='build/$PLATFORM/modellib', duplicate=0)
env.SConscript(dirs = ['doc'], build_dir='build/$PLATFORM/doc', duplicate=0)
env.SConscript(dirs = ['pyvisi/py_src'], build_dir='build/$PLATFORM/pyvisi', duplicate=0)
env.SConscript(dirs = ['pycad/py_src'], build_dir='build/$PLATFORM/pycad', duplicate=0)
env.SConscript(dirs = ['pythonMPI/src'], build_dir='build/$PLATFORM/pythonMPI', duplicate=0)

############ Targets to build and install libraries ############

target_init = env.Command(env['pyinstall']+'/__init__.py', None, Touch('$TARGET'))
env.Alias('target_init', [target_init])

# The headers have to be installed prior to build in order to satisfy #include <paso/Common.h>
env.Alias('build_esysUtils', ['target_install_esysUtils_headers', 'target_esysUtils_a'])
env.Alias('install_esysUtils', ['build_esysUtils', 'target_install_esysUtils_a'])

env.Alias('build_paso', ['target_install_paso_headers', 'target_paso_a'])
env.Alias('install_paso', ['build_paso', 'target_install_paso_a'])

env.Alias('build_escript', ['target_install_escript_headers', 'target_escript_so', 'target_escriptcpp_so'])
env.Alias('install_escript', ['build_escript', 'target_install_escript_so', 'target_install_escriptcpp_so', 'target_install_escript_py'])

env.Alias('build_finley', ['target_install_finley_headers', 'target_finley_so', 'target_finleycpp_so'])
env.Alias('install_finley', ['build_finley', 'target_install_finley_so', 'target_install_finleycpp_so', 'target_install_finley_py'])

# Now gather all the above into a couple easy targets: build_all and install_all
build_all_list = []
build_all_list += ['build_esysUtils']
build_all_list += ['build_paso']
build_all_list += ['build_escript']
build_all_list += ['build_finley']
if env['usempi']: build_all_list += ['target_pythonMPI_exe']
env.Alias('build_all', build_all_list)

install_all_list = []
install_all_list += ['target_init']
install_all_list += ['install_esysUtils']
install_all_list += ['install_paso']
install_all_list += ['install_escript']
install_all_list += ['install_finley']
install_all_list += ['target_install_pyvisi_py']
install_all_list += ['target_install_modellib_py']
install_all_list += ['target_install_pycad_py']
if env['usempi']: install_all_list += ['target_install_pythonMPI_exe']
env.Alias('install_all', install_all_list)

# Default target is install
env.Default('install_all')

############ Targets to build and run the test suite ###########

env.Alias('build_cppunittest', ['target_install_cppunittest_headers', 'target_cppunittest_a'])
env.Alias('install_cppunittest', ['build_cppunittest', 'target_install_cppunittest_a'])
env.Alias('run_tests', ['install_all', 'target_install_cppunittest_a'])
env.Alias('all_tests', ['install_all', 'target_install_cppunittest_a', 'run_tests', 'py_tests'])

############ Targets to build the documentation ################

env.Alias('docs', ['examples_tarfile', 'examples_zipfile', 'api_epydoc', 'api_doxygen', 'guide_pdf', 'guide_html'])

############ Remember what options were used in the build ######

if not IS_WINDOWS_PLATFORM:
  env.Execute("/bin/rm -f " + env['libinstall'] + "/Compiled.with.*")
  if env['usedebug']:		env.Execute("touch " + env['libinstall'] + "/Compiled.with.debug")
  if env['usempi']:		env.Execute("touch " + env['libinstall'] + "/Compiled.with.mpi")
  if env['omp_optim'] != '':	env.Execute("touch " + env['libinstall'] + "/Compiled.with.OpenMP")

