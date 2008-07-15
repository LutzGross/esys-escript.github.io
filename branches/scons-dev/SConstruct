#         Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php

# TODO use Install() to install python files instead of custom builder build_py?
# TODO How to modify CCFLAGS for only one file? Look for example with Program(...)
# TODO use USE_VTK in saveVTK to enable/disable this method
# TODO set CPPDEFINES and CXXDEFINES the same
# TODO configure for both 96.91 and 98.03

EnsureSConsVersion(0,96,91)
EnsurePythonVersion(2,3)

import sys, os, re, socket

# Add our extensions
if os.path.isdir('scons'): sys.path.append('scons')
import scons_extensions

# The string python2.4 or python2.5
python_version = 'python%s.%s' % (sys.version_info[0], sys.version_info[1])

# Which /usr/lib? Prefer /usr/lib64 if available for 64-bit systems
usr_lib = '/usr/lib'
if os.path.isdir('/usr/lib64'): usr_lib = '/usr/lib64'

# MS Windows support, many thanks to PH
IS_WINDOWS_PLATFORM = (os.name== "nt")

prefix = ARGUMENTS.get('prefix', Dir('#.').abspath)
print "Install prefix is: ", prefix

# Read configuration options from file scons/<hostname>_options.py
hostname = re.sub("[^0-9a-zA-Z]", "_", socket.gethostname().split('.')[0])
tmp = os.path.join("scons",hostname+"_options.py")
options_file = ARGUMENTS.get('options_file', tmp)
if not os.path.isfile(options_file): options_file = False

# Is setting options based on compiler name adequate instead of these?
# If scons/<hostname>_options.py not found try to specify defaults based on OS
#if not options_file:
#   if IS_WINDOWS_PLATFORM:
#      options_file = "scons/windows_mscv71_options.py"
#   else:
#      options_file = "scons/linux_gcc_eg_options.py"

# Load options file and command-line arguments
opts = Options(options_file, ARGUMENTS)

#==============================================================================================     
# Default options and options help text
# These are defaults and can be overridden using command line arguments or an options file.
# if the options_file or ARGUMENTS do not exist then the ones listed as default here are used
# DO NOT CHANGE THEM HERE
opts.AddOptions(
# Where to install esys stuff
  ('incinstall', 'where the esys headers will be installed',             Dir('#.').abspath+'/include'),
  ('libinstall', 'where the esys libraries will be installed',           os.path.join(prefix,"lib")),
  ('pyinstall', 'where the esys python modules will be installed',       os.path.join(prefix,"esys")),
# Compilation options
  BoolOption('dodebug', 'Do you want a debug build?', 'no'),
  ('options_file', "File of paths/options. Default: scons/<hostname>_options.py", options_file),
  ('cc_defines','C/C++ defines to use', None),

  ('cc_flags', 'C compiler flags to use', '-DEFAULT_1'),
  ('cp_flags', 'C++ compiler flags to use', '-DEFAULT_2'),
  ('cc_optim', 'C compiler optimization flags to use', '-DEFAULT_3'),
  ('cp_optim', 'C++ compiler optimization flags to use', '-DEFAULT_4'),
  ('cc_debug', 'C compiler debug flags to use', '-DEFAULT_5'),
  ('cp_debug', 'C++ compiler debug flags to use', '-DEFAULT_6'),
  ('omp_flags', 'OpenMP compiler flags to use (Release build)', '-DEFAULT_7'),
  ('omp_flags_debug', 'OpenMP compiler flags to use (Debug build)', '-DEFAULT_8'),
  ('sys_libs', 'System libraries to link with', []),
  ('ar_flags', 'Static library archiver flags to use', ''),
# Python
  PathOption('python_path', 'Path to Python includes', sys.prefix+'/include/'+python_version),
  PathOption('python_lib_path', 'Path to Python libs', sys.prefix+'/lib'),
  ('python_lib', 'Python libraries to link with', python_version),
  ('python_cmd', 'Python command', 'python'),
# Boost
  PathOption('boost_path', 'Path to Boost includes', "/usr/include"),
  PathOption('boost_lib_path', 'Path to Boost libs', usr_lib),
  ('boost_lib', 'Boost libraries to link with', ['boost_python']),
# netCDF
  ('useNetCDF', 'switch on/off the usage of netCDF', '1'),
  PathOption('netCDF_path', 'Path to netCDF includes', '/usr/include'),
  PathOption('netCDF_lib_path', 'Path to netCDF libs', usr_lib),
  ('netCDF_libs', 'netCDF C++ libraries to link with', ['netcdf_c++', 'netcdf']),
# MPI
  BoolOption('useMPI', 'Compile parallel version using MPI', 'no'),
  ('MPICH_IGNORE_CXX_SEEK', 'name of macro to ignore MPI settings of C++ SEEK macro (for MPICH)' , 'MPICH_IGNORE_CXX_SEEK'),
  PathOption('mpi_path', 'Path to MPI includes', '/usr/include'),
  ('mpi_run', 'mpirun name' , 'mpiexec -np 1'),
  PathOption('mpi_lib_path', 'Path to MPI libs (needs to be added to the LD_LIBRARY_PATH)', usr_lib),
  ('mpi_libs', 'MPI libraries to link with (needs to be shared!)', ['mpich' , 'pthread', 'rt']),
# ParMETIS
  BoolOption('useParMETIS', 'Compile parallel version using ParMETIS', 'no'),
  ('parmetis_path', 'Path to ParMETIS includes', '/usr/include'),
  ('parmetis_lib_path', 'Path to ParMETIS library', usr_lib),
  ('parmetis_libs', 'ParMETIS library to link with', []),
# PAPI
  PathOption('papi_path', 'Path to PAPI includes', None),
  PathOption('papi_lib_path', 'Path to PAPI libs', None),
  ('papi_libs', 'PAPI libraries to link with', []),
  ('papi_instrument_solver', 'use PAPI in Solver.c to instrument each iteration of the solver', False),
# MKL
  PathOption('mkl_path', 'Path to MKL includes', None),
  PathOption('mkl_lib_path', 'Path to MKL libs', None),
  ('mkl_libs', 'MKL libraries to link with', []),
# SCSL
  PathOption('scsl_path', 'Path to SCSL includes', None),
  PathOption('scsl_lib_path', 'Path to SCSL libs', None),
  ('scsl_libs', 'SCSL libraries to link with', []),
  ('scsl_libs_MPI', 'SCSL libraries to link with for MPI build', []),
# UMFPACK
  PathOption('ufc_path', 'Path to UFconfig includes', '/usr/include'),
  PathOption('umf_path', 'Path to UMFPACK includes', '/usr/include'),
  PathOption('umf_lib_path', 'Path to UMFPACK libs', usr_lib),
  ('umf_libs', 'UMFPACK libraries to link with', []),
# AMD (used by UMFPACK)
  PathOption('amd_path', 'Path to AMD includes', '/usr/include'),
  PathOption('amd_lib_path', 'Path to AMD libs', usr_lib),
  ('amd_libs', 'AMD libraries to link with', []),
# BLAS
  PathOption('blas_path', 'Path to BLAS includes', None),
  PathOption('blas_lib_path', 'Path to BLAS libs', None),
  ('blas_libs', 'BLAS libraries to link with', []),
)

# Specify which compilers to use (intelc uses regular expressions
# improperly and emits a warning about failing to find the compilers. 
# It can be safely ignored)

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

# Default compiler options (override allowed in hostname_options.py, but should not be necessary)
if env["CC"] == "icc":
  # Intel compilers
  cc_flags		= "-ansi -wd161 -w1 -vec-report0 -DBLOCKTIMER -DCORE_ID1"
  cp_flags		= "-ansi -wd161"
  cc_optim		= "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias"
  cp_optim		= "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias"
  cc_debug		= "-g -O0 -DDOASSERT -DDOPROF"
  cp_debug		= "-g -O0 -DDOASSERT -DDOPROF"
  omp_flags		= "-openmp -openmp_report0"
  omp_flags_debug	= "-openmp -openmp_report0"
elif env["CC"] == "gcc":
  # GNU C on any system
  cc_flags		= "-ffast-math -Wno-unknown-pragmas -ansi -DBLOCKTIMER -DCORE_ID1 -pedantic-errors -Wno-long-long"
  cp_flags		= "--no-warn -ansi"
  cc_optim		= "-O3"
  cp_optim		= "-O3"
  cc_debug		= "-g -O0 -DDOASSERT -DDOPROF"
  cp_debug		= "-g -O0 -DDOASSERT -DDOPROF"
  omp_flags		= ""
  omp_flags_debug	= ""
elif env["CC"] == "cl":
  # Microsoft Visual C on Windows
  cc_flags		= "/FD /EHsc /GR /wd4068 -D_USE_MATH_DEFINES -DDLL_NETCDF"
  cp_flags		= "/FD /EHsc /GR /wd4068"
  cc_optim		= "/O2 /Op /MT /W3"
  cp_optim		= "/O2 /Op /MT /W3"
  cc_debug		= "/Od /RTC1 /MTd /ZI"
  cp_debug		= "/Od /RTC1 /MTd /ZI"
  omp_flags		= ""
  omp_flags_debug	= ""

# If not specified in hostname_options.py then set them here
if env["cc_flags"]		== "-DEFAULT_1": env['cc_flags'] = cc_flags
if env["cp_flags"]		== "-DEFAULT_2": env['cp_flags'] = cp_flags
if env["cc_optim"]		== "-DEFAULT_3": env['cc_optim'] = cc_optim
if env["cp_optim"]		== "-DEFAULT_4": env['cp_optim'] = cp_optim
if env["cc_debug"]		== "-DEFAULT_5": env['cc_debug'] = cc_debug
if env["cp_debug"]		== "-DEFAULT_6": env['cp_debug'] = cp_debug
if env["omp_flags"]		== "-DEFAULT_7": env['omp_flags'] = omp_flags
if env["omp_flags_debug"]	== "-DEFAULT_8": env['omp_flags_debug'] = omp_flags_debug

# Get the global Subversion revision number for getVersion() method
try:
   global_revision = os.popen("svnversion -n .").read()
   global_revision = re.sub(":.*", "", global_revision)
   global_revision = re.sub("[^0-9]", "", global_revision)
except:
   global_revision="-1"
if global_revision == "": global_revision="-2"
env.Append(CPPDEFINES = ["SVN_VERSION="+global_revision])

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

# A few creature conveniences
libinstall = env['libinstall']
incinstall = env['incinstall']
pyinstall  = env['pyinstall']
sys_libs   = env['sys_libs']

# Configure for test suite: python run_testname.py
env.PrependENVPath('PYTHONPATH', prefix)
env.PrependENVPath('LD_LIBRARY_PATH', libinstall)

############ Set up paths for Configure() ######################

# Add gcc option -I<Escript>/trunk/include
env.Append(CPPPATH		= [Dir('include')])

# Add gcc option -L<Escript>/trunk/lib
env.Append(LIBPATH		= [Dir('lib')])

# Create a Configure() environment only for checking existence of
# libraries and headers.  Later we throw it away then build the real
# environment.
try:
  conf = Configure(env.Clone())	# scons-98.03
except AttributeError:
  conf = Configure(env.Copy())	# scons-96.92

# Must add all paths at first or Configure() can't cache results
conf.env.Append(CPPPATH	= [env['python_path']])
conf.env.Append(LIBS	= [env['python_lib']])
conf.env.Append(LIBPATH	= [env['python_lib_path']])
conf.env.Append(CPPPATH	= [env['boost_path']])
conf.env.Append(LIBS	= [env['boost_lib']])
conf.env.Append(LIBPATH	= [env['boost_lib_path']])
conf.env.Append(CPPPATH	= [env['netCDF_path']])
conf.env.Append(LIBS	= [env['netCDF_libs']])
conf.env.Append(LIBPATH	= [env['netCDF_lib_path']])
conf.env.Append(CPPPATH	= [env['mpi_path']])
conf.env.Append(LIBS	= [env['mpi_libs']])
conf.env.Append(LIBPATH	= [env['mpi_lib_path']])
conf.env.Append(CPPPATH	= [env['parmetis_path']])
conf.env.Append(LIBS	= [env['parmetis_libs']])
conf.env.Append(LIBPATH	= [env['parmetis_lib_path']])

############ numarray (required) ###############################

try: from numarray import identity
except ImportError: sys.exit(1)

############ python (required) #################################

if not conf.CheckCHeader('Python.h'): sys.exit(1)
if not conf.CheckFunc('Py_Main'): sys.exit(1)

############ boost (required) ##################################

if not conf.CheckCXXHeader('boost/python.hpp'): sys.exit(1)
if not conf.CheckFunc('PyObject_SetAttr'): sys.exit(1)

############ NetCDF (optional) #################################

env['useNetCDF'] = 1

if env['useNetCDF'] and not conf.CheckCHeader('netcdf.h'): env['useNetCDF'] = 0
if env['useNetCDF'] and not conf.CheckFunc('nc_open'): env['useNetCDF'] = 0

############ VTK (optional) ####################################

try:
  import vtk
  useVTK = 1
except ImportError:
  useVTK = 0

############ MPI (optional) ####################################

if env['useMPI'] and not conf.CheckCHeader('mpi.h'): env['useMPI'] = 0
if env['useMPI'] and not conf.CheckFunc('MPI_Init'): env['useMPI'] = 0

############ ParMETIS (optional) ###############################

env['useParMETIS'] = 1
if not env['useMPI']: env['useParMETIS'] = 0
if env['useParMETIS'] and not conf.CheckCHeader('parmetis.h'): env['useParMETIS'] = 0
if env['useParMETIS'] and not conf.CheckFunc('ParMETIS_V3_PartGeomKway'): env['useParMETIS'] = 0

############ Configure finished, now modify environment ########

conf.Finish()

# Enable debug
if env['dodebug']:
  env.Append(CCFLAGS		= env['cc_debug'])
  env.Append(CCFLAGS		= env['omp_flags_debug'])
  env.Append(CCFLAGS		= env['cc_flags'])
  env.Append(CPPDEFINES		= ['BOUNDS_CHECK'])
else:
  env.Append(CCFLAGS		= env['cc_optim'])
  env.Append(CCFLAGS		= env['omp_flags'])
  env.Append(CCFLAGS		= env['cc_flags'])

# Python
env.Append(CPPPATH		= [env['python_path']])
env.Append(LIBPATH		= [env['python_lib_path']])
env.Append(LIBS			= [env['python_lib']])

# Boost
env.Append(CPPPATH		= [env['boost_path']])
env.Append(LIBPATH		= [env['boost_lib_path']])
env.Append(LIBS			= [env['boost_lib']])

# NetCDF
if env['useNetCDF']:
  env.Append(CPPPATH		= [env['netCDF_path']])
  env.Append(LIBPATH		= [env['netCDF_lib_path']])
  env.Append(LIBS		= [env['netCDF_libs']])
  env.Append(CPPDEFINES		= ['USE_NETCDF'])

# VTK
if useVTK:
  env.Append(CPPDEFINES		= ['USE_VTK'])

# MS Windows
if IS_WINDOWS_PLATFORM:
  env.PrependENVPath('PATH',	[env['boost_lib_path']])
  env.PrependENVPath('PATH',	libinstall)
  if env['useNetCDF']:
    env.PrependENVPath('PATH',	[env['netCDF_lib_path']])

# Create a modified environment for MPI programs
try:
  env_mpi = env.Clone()	# scons-98.03
except AttributeError:
  env_mpi = env.Copy()	# scons-96.92

# MPI
if env_mpi['useMPI']:
  env_mpi.Append(CPPPATH	= [env['mpi_path']])
  env_mpi.Append(LIBPATH	= [env['mpi_lib_path']])
  env_mpi.Append(LIBS		= [env['mpi_libs']])
  env_mpi.Append(CPPDEFINES	= ['PASO_MPI', 'MPI_NO_CPPBIND', env_mpi['MPICH_IGNORE_CXX_SEEK']])

# ParMETIS
if env_mpi['useParMETIS']:
  env_mpi.Append(CPPPATH	= [env['parmetis_path']])
  env_mpi.Append(LIBPATH	= [env['parmetis_lib_path']])
  env_mpi.Append(LIBS		= [env['parmetis_libs']])
  env_mpi.Append(CPPDEFINES	= ['USE_PARMETIS'])

env.Append(ARFLAGS = env['ar_flags'])

############ End of modify environment #########################

#==========================================================================
#
#    Add some custom builders
#
py_builder = Builder(action = scons_extensions.build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed',
                              src_suffix=env['PROGSUFFIX'], single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = scons_extensions.runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

# ============= Remember what options were used in the compile =====================================
if not IS_WINDOWS_PLATFORM:
  env.Execute("/bin/rm -f " + libinstall + "/Compiled.with.*")
  if env['dodebug']:		env.Execute("touch " + libinstall + "/Compiled.with.debug")
  if env['useMPI']:		env.Execute("touch " + libinstall + "/Compiled.with.mpi")
  if env['omp_flags'] != '':	env.Execute("touch " + libinstall + "/Compiled.with.OpenMP")

Export(["env", "env_mpi", "incinstall", "libinstall", "pyinstall", "sys_libs", "prefix" ])

env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
# env.SConscript(dirs = ['bruce/src'], build_dir='build/$PLATFORM/bruce', duplicate=0)
env.SConscript(dirs = ['escript/src'], build_dir='build/$PLATFORM/escript', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['finley/src'], build_dir='build/$PLATFORM/finley', duplicate=0)
env.SConscript(dirs = ['modellib/py_src'], build_dir='build/$PLATFORM/modellib', duplicate=0)
env.SConscript(dirs = ['doc'], build_dir='build/$PLATFORM/doc', duplicate=0)
env.SConscript(dirs = ['pyvisi/py_src'], build_dir='build/$PLATFORM/pyvisi', duplicate=0)
env.SConscript(dirs = ['pycad/py_src'], build_dir='build/$PLATFORM/pycad', duplicate=0)
env.SConscript(dirs = ['pythonMPI/src'], build_dir='build/$PLATFORM/pythonMPI', duplicate=0)

############ Targets to build and install libraries ############

target_init = env.Command(pyinstall+'/__init__.py', None, Touch('$TARGET'))
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
# build_all_list += ['build_bruce']

install_all_list = []
install_all_list += ['target_init']
install_all_list += ['install_esysUtils']
install_all_list += ['install_paso']
install_all_list += ['install_escript']
install_all_list += ['install_finley']
install_all_list += ['target_install_pyvisi_py']
install_all_list += ['target_install_modellib_py']
install_all_list += ['target_install_pycad_py']
# install_all_list += ['target_install_bruce']
# install_all_list += ['target_install_bruce_py']

env.Alias('build_all', build_all_list)
env.Alias('install_all', install_all_list)

env.Default('install_all')

############ Targets to build and run the test suite ###########

env.Alias('build_cppunittest', ['target_install_cppunittest_headers', 'target_cppunittest_a'])
env.Alias('install_cppunittest', ['build_cppunittest', 'target_install_cppunittest_a'])
env.Alias('run_tests', ['install_all', 'target_install_cppunittest_a'])
env.Alias('all_tests', ['install_all', 'target_install_cppunittest_a', 'run_tests', 'py_tests'])

############ Targets to build the documentation ################

env.Alias('docs', ['examples_tarfile', 'examples_zipfile', 'api_epydoc', 'api_doxygen', 'guide_pdf', 'guide_html'])

