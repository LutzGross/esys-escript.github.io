
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################


EnsureSConsVersion(0,96,91)
EnsurePythonVersion(2,3)

import sys, os, re, socket, platform, stat

# Add our extensions
if os.path.isdir('scons'): sys.path.append('scons')
import scons_extensions

# Use /usr/lib64 if available, else /usr/lib
usr_lib = '/usr/lib'
if os.path.isfile('/usr/lib64/libc.so'): usr_lib = '/usr/lib64'

# The string python2.4 or python2.5
python_version = 'python%s.%s' % (sys.version_info[0], sys.version_info[1])

# MS Windows support, many thanks to PH
IS_WINDOWS_PLATFORM = (os.name== "nt")

prefix = ARGUMENTS.get('prefix', Dir('#.').abspath)

#Determine where to read options from use:
#1. command line
#2. scons/<hostname>_options.py
#3. name as part of a cluster
options_file=ARGUMENTS.get('options_file', None)
if not options_file:
  hostname = re.sub("[^0-9a-zA-Z]", "_", socket.gethostname().split('.')[0])
  options_file = os.path.join("scons",hostname+"_options.py")
  #If there is no options file with that name see if there is a substitute
  if not os.path.isfile(options_file):
    tmp = scons_extensions.effectiveName(hostname)
    options_file = os.path.join("scons",tmp+"_options.py")

if not os.path.isfile(options_file):
  print "Options file not found (expected '%s')" % options_file
  options_file = False
else:
  print "Options file is", options_file

# Load options file and command-line arguments
opts = Options(options_file, ARGUMENTS)

############ Load build options ################################

opts.AddOptions(
# Where to install esys stuff
  ('prefix', 'where everything will be installed',                       Dir('#.').abspath),
  ('incinstall', 'where the esys headers will be installed',             os.path.join(Dir('#.').abspath,'include')),
  ('bininstall', 'where the esys binaries will be installed',            os.path.join(prefix,'bin')),
  ('libinstall', 'where the esys libraries will be installed',           os.path.join(prefix,'lib')),
  ('pyinstall', 'where the esys python modules will be installed',       os.path.join(prefix,'esys')),
# Compilation options
  BoolOption('dodebug', 'For backwards compatibility', 'no'),
  BoolOption('usedebug', 'Do you want a debug build?', 'no'),
  BoolOption('usevtk', 'Do you want to use VTK?', 'yes'),
  ('options_file', 'File of paths/options. Default: scons/<hostname>_options.py', options_file),
  ('win_cc_name', 'windows C compiler name if needed', 'msvc'),
  # The strings -DDEFAULT_ get replaced by scons/<hostname>_options.py or by defaults below
  ('cc_flags', 'C compiler flags to use', '-DEFAULT_1'),
  ('cc_optim', 'C compiler optimization flags to use', '-DEFAULT_2'),
  ('cc_debug', 'C compiler debug flags to use', '-DEFAULT_3'),
  ('omp_optim', 'OpenMP compiler flags to use (Release build)', '-DEFAULT_4'),
  ('omp_debug', 'OpenMP compiler flags to use (Debug build)', '-DEFAULT_5'),
  ('omp_libs', 'OpenMP compiler libraries to link with', '-DEFAULT_6'),
  ('cc_extra', 'Extra C/C++ flags', ''),
  ('ld_extra', 'Extra linker flags', ''),
  ('sys_libs', 'System libraries to link with', []),
  ('ar_flags', 'Static library archiver flags to use', ''),
  BoolOption('useopenmp', 'Compile parallel version using OpenMP', 'no'),
  BoolOption('usepedantic', 'Compile with -pedantic if using gcc', 'no'),
  BoolOption('usewarnings','Compile with warnings as errors if using gcc','yes'),
  ('forcelazy','for testing use only - set the default value for autolazy','leave_alone'),
# Python
  ('python_path', 'Path to Python includes', '/usr/include/'+python_version),
  ('python_lib_path', 'Path to Python libs', usr_lib),
  ('python_libs', 'Python libraries to link with', [python_version]),
  ('python_cmd', 'Python command', 'python'),
# Boost
  ('boost_path', 'Path to Boost includes', '/usr/include'),
  ('boost_lib_path', 'Path to Boost libs', usr_lib),
  ('boost_libs', 'Boost libraries to link with', ['boost_python']),
# NetCDF
  BoolOption('usenetcdf', 'switch on/off the usage of netCDF', 'yes'),
  ('netCDF_path', 'Path to netCDF includes', '/usr/include'),
  ('netCDF_lib_path', 'Path to netCDF libs', usr_lib),
  ('netCDF_libs', 'netCDF C++ libraries to link with', ['netcdf_c++', 'netcdf']),
# MPI
  BoolOption('useMPI', 'For backwards compatibility', 'no'),
  BoolOption('usempi', 'Compile parallel version using MPI', 'no'),
  ('MPICH_IGNORE_CXX_SEEK', 'name of macro to ignore MPI settings of C++ SEEK macro (for MPICH)' , 'MPICH_IGNORE_CXX_SEEK'),
  ('mpi_path', 'Path to MPI includes', '/usr/include'),
  ('mpi_run', 'mpirun name' , 'mpiexec -np 1'),
  ('mpi_lib_path', 'Path to MPI libs (needs to be added to the LD_LIBRARY_PATH)', usr_lib),
  ('mpi_libs', 'MPI libraries to link with (needs to be shared!)', ['mpich' , 'pthread', 'rt']),
  ('mpi_flavour','Type of MPI execution environment','none'), 
# ParMETIS
  BoolOption('useparmetis', 'Compile parallel version using ParMETIS', 'yes'),
  ('parmetis_path', 'Path to ParMETIS includes', '/usr/include'),
  ('parmetis_lib_path', 'Path to ParMETIS library', usr_lib),
  ('parmetis_libs', 'ParMETIS library to link with', ['parmetis', 'metis']),
# PAPI
  BoolOption('usepapi', 'switch on/off the usage of PAPI', 'no'),
  ('papi_path', 'Path to PAPI includes', '/usr/include'),
  ('papi_lib_path', 'Path to PAPI libs', usr_lib),
  ('papi_libs', 'PAPI libraries to link with', ['papi']),
  BoolOption('papi_instrument_solver', 'use PAPI in Solver.c to instrument each iteration of the solver', False),
# MKL
  BoolOption('usemkl', 'switch on/off the usage of MKL', 'no'),
  ('mkl_path', 'Path to MKL includes', '/sw/sdev/cmkl/10.0.2.18/include'),
  ('mkl_lib_path', 'Path to MKL libs', '/sw/sdev/cmkl/10.0.2.18/lib/em64t'),
  ('mkl_libs', 'MKL libraries to link with', ['mkl_solver', 'mkl_em64t', 'guide', 'pthread']),
# UMFPACK
  BoolOption('useumfpack', 'switch on/off the usage of UMFPACK', 'no'),
  ('ufc_path', 'Path to UFconfig includes', '/usr/include/suitesparse'),
  ('umf_path', 'Path to UMFPACK includes', '/usr/include/suitesparse'),
  ('umf_lib_path', 'Path to UMFPACK libs', usr_lib),
  ('umf_libs', 'UMFPACK libraries to link with', ['umfpack']),
# Silo
  BoolOption('usesilo', 'switch on/off the usage of Silo', 'yes'),
  ('silo_path', 'Path to Silo includes', '/usr/include'),
  ('silo_lib_path', 'Path to Silo libs', usr_lib),
  ('silo_libs', 'Silo libraries to link with', ['siloh5', 'hdf5']),
# AMD (used by UMFPACK)
  ('amd_path', 'Path to AMD includes', '/usr/include/suitesparse'),
  ('amd_lib_path', 'Path to AMD libs', usr_lib),
  ('amd_libs', 'AMD libraries to link with', ['amd']),
# BLAS (used by UMFPACK)
  ('blas_path', 'Path to BLAS includes', '/usr/include/suitesparse'),
  ('blas_lib_path', 'Path to BLAS libs', usr_lib),
  ('blas_libs', 'BLAS libraries to link with', ['blas']),
# An option for specifying the compiler tools set (see windows branch).
  ('tools_names', 'allow control over the tools in the env setup', ['intelc']),
# finer control over library building, intel aggressive global optimisation
# works with dynamic libraries on windows.
  ('share_esysUtils', 'control static or dynamic esysUtils lib', False),
  ('share_paso', 'control static or dynamic paso lib', False)
)

############ Specify which compilers to use ####################

# intelc uses regular expressions improperly and emits a warning about
# failing to find the compilers.  This warning can be safely ignored.

if IS_WINDOWS_PLATFORM:
      env = Environment(options = opts)
      env = Environment(tools = ['default'] + env['tools_names'],
                        options = opts)
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

sysheaderopt = ""		# how do we indicate that a header is a system header. Use "" for no action.

if env["CC"] == "icc":
  # Intel compilers
  cc_flags		= "-fPIC -ansi -wd161 -w1 -vec-report0 -DBLOCKTIMER -DCORE_ID1"
  cc_optim		= "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias"
  cc_debug		= "-g -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK"
  omp_optim		= "-openmp -openmp_report0"
  omp_debug		= "-openmp -openmp_report0"
  omp_libs		= ['guide', 'pthread']
  pedantic		= ""
  fatalwarning		= ""		# Switch to turn warnings into errors
  sysheaderopt		= ""
elif env["CC"] == "gcc":
  # GNU C on any system
  cc_flags		= "-pedantic -Wall -fPIC -ansi -ffast-math -Wno-unknown-pragmas -DBLOCKTIMER  -Wno-sign-compare -Wno-system-headers -Wno-long-long -Wno-strict-aliasing"
#the long long warning occurs on the Mac
  cc_optim		= "-O3"
  cc_debug		= "-g -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK"
  omp_optim		= "-fopenmp"
  omp_debug		= "-fopenmp"
  omp_libs		= ['gomp']
  pedantic		= "-pedantic-errors -Wno-long-long"
  fatalwarning		= "-Werror"
  sysheaderopt		= "-isystem "
elif env["CC"] == "cl":
  # Microsoft Visual C on Windows
  cc_flags		= "/FD /EHsc /GR /wd4068 -D_USE_MATH_DEFINES -DDLL_NETCDF"
  cc_optim		= "/O2 /Op /MT /W3"
  cc_debug		= "/Od /RTC1 /MTd /ZI -DBOUNDS_CHECK"
  omp_optim		= ""
  omp_debug		= ""
  omp_libs		= []
  pedantic		= ""
  fatalwarning		= ""
  sysheaderopt		= ""
elif env["CC"] == "icl":
  # intel C on Windows, see windows_intelc_options.py for a start
  pedantic		= ""
  fatalwarning		= ""
  sysheaderopt		= ""


# If not specified in hostname_options.py then set them here
if env["cc_flags"]	== "-DEFAULT_1": env['cc_flags'] = cc_flags
if env["cc_optim"]	== "-DEFAULT_2": env['cc_optim'] = cc_optim
if env["cc_debug"]	== "-DEFAULT_3": env['cc_debug'] = cc_debug
if env["omp_optim"]	== "-DEFAULT_4": env['omp_optim'] = omp_optim
if env["omp_debug"]	== "-DEFAULT_5": env['omp_debug'] = omp_debug
if env["omp_libs"]	== "-DEFAULT_6": env['omp_libs'] = omp_libs

#set up the autolazy values
if env['forcelazy']    != "leave_alone":
  if env['forcelazy'] == 'on':
	env.Append(CPPDEFINES='FAUTOLAZYON')
  else:
     if env['forcelazy'] == 'off':
	env.Append(CPPDEFINES='FAUTOLAZYOFF')

# OpenMP is disabled if useopenmp=no or both variables omp_optim and omp_debug are empty
if not env["useopenmp"]:
  env['omp_optim'] = ""
  env['omp_debug'] = ""
  env['omp_libs'] = []

if env['omp_optim'] == "" and env['omp_debug'] == "": env["useopenmp"] = 0

############ Copy environment variables into scons env #########

try: env['ENV']['OMP_NUM_THREADS'] = os.environ['OMP_NUM_THREADS']
except KeyError: env['ENV']['OMP_NUM_THREADS'] = 1

try: env['ENV']['ESCRIPT_NUM_THREADS'] = os.environ['ESCRIPT_NUM_THREADS']
except KeyError: pass

try: env['ENV']['ESCRIPT_NUM_PROCS'] = os.environ['ESCRIPT_NUM_PROCS']
except KeyError: env['ENV']['ESCRIPT_NUM_PROCS']=1

try: env['ENV']['ESCRIPT_NUM_NODES'] = os.environ['ESCRIPT_NUM_NODES']
except KeyError: env['ENV']['ESCRIPT_NUM_NODES']=1

try: env['ENV']['ESCRIPT_HOSTFILE'] = os.environ['ESCRIPT_HOSTFILE']
except KeyError: pass

try: env['ENV']['PATH'] = os.environ['PATH']
except KeyError: pass

try: env['ENV']['PYTHONPATH'] = os.environ['PYTHONPATH']
except KeyError: pass

try: env['ENV']['C_INCLUDE_PATH'] = os.environ['C_INCLUDE_PATH']
except KeyError: pass

try: env['ENV']['CPLUS_INCLUDE_PATH'] = os.environ['CPLUS_INCLUDE_PATH']
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

env['ENV']['ESCRIPT_ROOT'] = prefix

############ Set up paths for Configure() ######################

# Make a copy of an environment
# Use env.Clone if available, but fall back on env.Copy for older version of scons
def clone_env(env):
  if 'Clone' in dir(env): return env.Clone()	# scons-0.98
  else:                   return env.Copy()	# scons-0.96

# Add cc option -I<Escript>/trunk/include
env.Append(CPPPATH		= [Dir('include')])

# Add cc option -L<Escript>/trunk/lib
env.Append(LIBPATH		= [Dir(env['libinstall'])])

if env['cc_extra'] != '': env.Append(CCFLAGS = env['cc_extra'])
if env['ld_extra'] != '': env.Append(LINKFLAGS = env['ld_extra'])

if env['usepedantic']: env.Append(CCFLAGS = pedantic)

# MS Windows
if IS_WINDOWS_PLATFORM:
  env.AppendENVPath('PATH',	[env['boost_lib_path']])
  env.AppendENVPath('PATH',	[env['libinstall']])
  if not env['share_esysUtils'] :
    env.Append(CPPDEFINES = ['ESYSUTILS_STATIC_LIB'])
  if not env['share_paso'] :
    env.Append(CPPDEFINES = ['PASO_STATIC_LIB'])

  if env['usenetcdf']:
    env.AppendENVPath('PATH',	[env['netCDF_lib_path']])

env.Append(ARFLAGS = env['ar_flags'])

# Get the global Subversion revision number for getVersion() method
try:
   global_revision = os.popen("svnversion -n .").read()
   global_revision = re.sub(":.*", "", global_revision)
   global_revision = re.sub("[^0-9]", "", global_revision)
except:
   global_revision="-1"
if global_revision == "": global_revision="-2"
env.Append(CPPDEFINES = ["SVN_VERSION="+global_revision])

############ numarray (required) ###############################

try:
  from numarray import identity
except ImportError:
  print "Cannot import numarray, you need to set your PYTHONPATH"
  sys.exit(1)

############ C compiler (required) #############################

# Create a Configure() environment for checking existence of required libraries and headers
conf = Configure(clone_env(env))

# Test that the compiler is working
if not conf.CheckFunc('printf'):
   print "Cannot run C compiler '%s' (or libc is missing)" % (env['CC'])
   sys.exit(1)

if conf.CheckFunc('gethostname'):
  conf.env.Append(CPPDEFINES = ['HAVE_GETHOSTNAME'])

############ python libraries (required) #######################


if not sysheaderopt =="":
  conf.env.Append(CCFLAGS=sysheaderopt+env['python_path'])
else:
  conf.env.AppendUnique(CPPPATH		= [env['python_path']])

conf.env.AppendUnique(LIBPATH		= [env['python_lib_path']])
conf.env.AppendUnique(LIBS		= [env['python_libs']])

conf.env.PrependENVPath('LD_LIBRARY_PATH', env['python_lib_path'])	# The wrapper script needs to find these libs
conf.env.PrependENVPath('PYTHONPATH', prefix)
conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if not conf.CheckCHeader('Python.h'):
  print "Cannot find python include files (tried 'Python.h' in directory %s)" % (env['python_path'])
  sys.exit(1)
if not conf.CheckFunc('Py_Exit'):
  print "Cannot find python library method Py_Main (tried lib %s in directory %s)" % (env['python_libs'], env['python_lib_path'])
  sys.exit(1)

############ boost (required) ##################################

if not sysheaderopt =="":
# This is required because we can't -isystem /usr/system because it breaks std includes
  if os.path.normpath(env['boost_path']) =="/usr/include":
    conf.env.Append(CCFLAGS=sysheaderopt+os.path.join(env['boost_path'],'boost'))
  else:
    conf.env.Append(CCFLAGS=sysheaderopt+env['boost_path'])
else:
  conf.env.AppendUnique(CPPPATH		= [env['boost_path']])

conf.env.AppendUnique(LIBPATH		= [env['boost_lib_path']])
conf.env.AppendUnique(LIBS		= [env['boost_libs']])

conf.env.PrependENVPath('LD_LIBRARY_PATH', env['boost_lib_path'])	# The wrapper script needs to find these libs
#ensure that our path entries remain at the front
conf.env.PrependENVPath('PYTHONPATH', prefix)
conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if not conf.CheckCXXHeader('boost/python.hpp'):
  print "Cannot find boost include files (tried boost/python.hpp in directory %s)" % (env['boost_path'])
  sys.exit(1)

if not conf.CheckFunc('PyObject_SetAttr'):
  print "Cannot find boost library method PyObject_SetAttr (tried method PyObject_SetAttr in library %s in directory %s)" % (env['boost_libs'], env['boost_lib_path'])
  sys.exit(1)

# Commit changes to environment
env = conf.Finish()

############ VTK (optional) ####################################

if env['usevtk']:
  try:
    import vtk
    env['usevtk'] = 1
  except ImportError:
    env['usevtk'] = 0

# Add VTK to environment env if it was found
if env['usevtk']:
  env.Append(CPPDEFINES = ['USE_VTK'])

############ NetCDF (optional) #################################

conf = Configure(clone_env(env))

if env['usenetcdf']:
  conf.env.AppendUnique(CPPPATH	= [env['netCDF_path']])
  conf.env.AppendUnique(LIBPATH	= [env['netCDF_lib_path']])
  conf.env.AppendUnique(LIBS	= [env['netCDF_libs']])
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['netCDF_lib_path'])	# The wrapper script needs to find these libs
  #ensure that our path entries remain at the front
  conf.env.PrependENVPath('PYTHONPATH', prefix)
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if env['usenetcdf'] and not conf.CheckCHeader('netcdf.h'): env['usenetcdf'] = 0
if env['usenetcdf'] and not conf.CheckFunc('nc_open'): env['usenetcdf'] = 0

# Add NetCDF to environment env if it was found
if env['usenetcdf']:
  env = conf.Finish()
  env.Append(CPPDEFINES = ['USE_NETCDF'])
else:
  conf.Finish()

############ PAPI (optional) ###################################

# Start a new configure environment that reflects what we've already found
conf = Configure(clone_env(env))

if env['usepapi']:
  conf.env.AppendUnique(CPPPATH	= [env['papi_path']])
  conf.env.AppendUnique(LIBPATH	= [env['papi_lib_path']])
  conf.env.AppendUnique(LIBS	= [env['papi_libs']])
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['papi_lib_path'])	# The wrapper script needs to find these libs
  #ensure that our path entries remain at the front
  conf.env.PrependENVPath('PYTHONPATH', prefix)
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if env['usepapi'] and not conf.CheckCHeader('papi.h'): env['usepapi'] = 0
if env['usepapi'] and not conf.CheckFunc('PAPI_start_counters'): env['usepapi'] = 0

# Add PAPI to environment env if it was found
if env['usepapi']:
  env = conf.Finish()
  env.Append(CPPDEFINES = ['BLOCKPAPI'])
else:
  conf.Finish()

############ MKL (optional) ####################################

# Start a new configure environment that reflects what we've already found
conf = Configure(clone_env(env))

if env['usemkl']:
  conf.env.AppendUnique(CPPPATH	= [env['mkl_path']])
  conf.env.AppendUnique(LIBPATH	= [env['mkl_lib_path']])
  conf.env.AppendUnique(LIBS	= [env['mkl_libs']])
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['mkl_lib_path'])	# The wrapper script needs to find these libs
  #ensure that our path entries remain at the front
  conf.env.PrependENVPath('PYTHONPATH', prefix)
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if env['usemkl'] and not conf.CheckCHeader('mkl_solver.h'): env['usemkl'] = 0
if env['usemkl'] and not conf.CheckFunc('pardiso'): env['usemkl'] = 0

# Add MKL to environment env if it was found
if env['usemkl']:
  env = conf.Finish()
  env.Append(CPPDEFINES = ['MKL'])
else:
  conf.Finish()

############ UMFPACK (optional) ################################

# Start a new configure environment that reflects what we've already found
conf = Configure(clone_env(env))

if env['useumfpack']:
  conf.env.AppendUnique(CPPPATH	= [env['ufc_path']])
  conf.env.AppendUnique(CPPPATH	= [env['umf_path']])
  conf.env.AppendUnique(LIBPATH	= [env['umf_lib_path']])
  conf.env.AppendUnique(LIBS	= [env['umf_libs']])
  conf.env.AppendUnique(CPPPATH	= [env['amd_path']])
  conf.env.AppendUnique(LIBPATH	= [env['amd_lib_path']])
  conf.env.AppendUnique(LIBS	= [env['amd_libs']])
  conf.env.AppendUnique(CPPPATH	= [env['blas_path']])
  conf.env.AppendUnique(LIBPATH	= [env['blas_lib_path']])
  conf.env.AppendUnique(LIBS	= [env['blas_libs']])
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['umf_lib_path'])	# The wrapper script needs to find these libs
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['amd_lib_path'])	# The wrapper script needs to find these libs
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['blas_lib_path'])	# The wrapper script needs to find these libs
  #ensure that our path entries remain at the front
  conf.env.PrependENVPath('PYTHONPATH', prefix)
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if env['useumfpack'] and not conf.CheckFunc('umfpack_di_symbolic'): env['useumfpack'] = 0
if env['useumfpack'] and not conf.CheckCHeader('umfpack.h'): env['useumfpack'] = 0
# if env['useumfpack'] and not conf.CheckFunc('daxpy'): env['useumfpack'] = 0 # this does not work on shake73?

# Add UMFPACK to environment env if it was found
if env['useumfpack']:
  env = conf.Finish()
  env.Append(CPPDEFINES = ['UMFPACK'])
else:
  conf.Finish()

############ Silo (optional) ###################################

if env['usesilo']:
  conf = Configure(clone_env(env))
  conf.env.AppendUnique(CPPPATH = [env['silo_path']])
  conf.env.AppendUnique(LIBPATH = [env['silo_lib_path']])
  conf.env.AppendUnique(LIBS = [env['silo_libs']])
  if not conf.CheckCHeader('silo.h'): env['usesilo'] = 0
  if not conf.CheckFunc('DBMkDir'): env['usesilo'] = 0
  conf.Finish()

# Add the path to Silo to environment env if it was found.
# Note that we do not add the libs since they are only needed for the
# escriptreader library and tools.
if env['usesilo']:
  env.AppendUnique(CPPPATH = [env['silo_path']])
  env.AppendUnique(LIBPATH = [env['silo_lib_path']])
  env.Append(CPPDEFINES = ['HAVE_SILO'])

############ Add the compiler flags ############################

# Enable debug by choosing either cc_debug or cc_optim
if env['usedebug']:
  env.Append(CCFLAGS		= env['cc_debug'])
  env.Append(CCFLAGS		= env['omp_debug'])
else:
  env.Append(CCFLAGS		= env['cc_optim'])
  env.Append(CCFLAGS		= env['omp_optim'])

# Always use cc_flags
env.Append(CCFLAGS		= env['cc_flags'])
env.Append(LIBS			= [env['omp_libs']])

############ Add some custom builders ##########################

py_builder = Builder(action = scons_extensions.build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed', src_suffix=env['PROGSUFFIX'], single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = scons_extensions.runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

epstopdfbuilder = Builder(action = scons_extensions.eps2pdf, suffix=".pdf", src_suffix=".eps", single_source=True)
env.Append(BUILDERS = {'EpsToPDF' : epstopdfbuilder});

############ MPI (optional) ####################################
if not env['usempi']: env['mpi_flavour']='none'

# Create a modified environment for MPI programs (identical to env if usempi=no)
env_mpi = clone_env(env)

# Start a new configure environment that reflects what we've already found
conf = Configure(clone_env(env_mpi))

if env_mpi['usempi']:
  VALID_MPIs=[ "MPT", "MPICH", "MPICH2", "OPENMPI", "INTELMPI" ]
  if not env_mpi['mpi_flavour'] in VALID_MPIs: 
      raise ValueError,"MPI is enabled but mpi_flavour = %s is not a valid key from %s."%( env_mpi['mpi_flavour'],VALID_MPIs)
  conf.env.AppendUnique(CPPPATH	= [env_mpi['mpi_path']])
  conf.env.AppendUnique(LIBPATH	= [env_mpi['mpi_lib_path']])
  conf.env.AppendUnique(LIBS	= [env_mpi['mpi_libs']])
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['mpi_lib_path'])	# The wrapper script needs to find these libs
  #ensure that our path entries remain at the front
  conf.env.PrependENVPath('PYTHONPATH', prefix)
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if env_mpi['usempi'] and not conf.CheckCHeader('mpi.h'): env_mpi['usempi'] = 0
# if env_mpi['usempi'] and not conf.CheckFunc('MPI_Init'): env_mpi['usempi'] = 0

# Add MPI to environment env_mpi if it was found
if env_mpi['usempi']:
  env_mpi = conf.Finish()
  env_mpi.Append(CPPDEFINES = ['PASO_MPI', 'MPI_NO_CPPBIND', env_mpi['MPICH_IGNORE_CXX_SEEK']])
else:
  conf.Finish()

env['usempi'] = env_mpi['usempi']


############ ParMETIS (optional) ###############################

# Start a new configure environment that reflects what we've already found
conf = Configure(clone_env(env_mpi))

if not env_mpi['usempi']: env_mpi['useparmetis'] = 0

if env_mpi['useparmetis']:
  conf.env.AppendUnique(CPPPATH	= [env_mpi['parmetis_path']])
  conf.env.AppendUnique(LIBPATH	= [env_mpi['parmetis_lib_path']])
  conf.env.AppendUnique(LIBS	= [env_mpi['parmetis_libs']])
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['parmetis_lib_path'])	# The wrapper script needs to find these libs
  #ensure that our path entries remain at the front
  conf.env.PrependENVPath('PYTHONPATH', prefix)
  conf.env.PrependENVPath('LD_LIBRARY_PATH', env['libinstall'])

if env_mpi['useparmetis'] and not conf.CheckCHeader('parmetis.h'): env_mpi['useparmetis'] = 0
if env_mpi['useparmetis'] and not conf.CheckFunc('ParMETIS_V3_PartGeomKway'): env_mpi['useparmetis'] = 0

# Add ParMETIS to environment env_mpi if it was found
if env_mpi['useparmetis']:
  env_mpi = conf.Finish()
  env_mpi.Append(CPPDEFINES = ['USE_PARMETIS'])
else:
  conf.Finish()

env['useparmetis'] = env_mpi['useparmetis']

############ Now we switch on Warnings as errors ###############

#this needs to be done after configuration because the scons test files have warnings in them

if ((fatalwarning != "") and (env['usewarnings'])):
  env.Append(CCFLAGS		= fatalwarning)
  env_mpi.Append(CCFLAGS		= fatalwarning)

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
if env['usesilo']: print "	Using Silo"
else: print "	Not using Silo"
if env['useopenmp']: print "	Using OpenMP"
else: print "	Not using OpenMP"
if env['usempi']: print "	Using MPI (flavour = %s)"%env['mpi_flavour']
else: print "	Not using MPI"
if env['useparmetis']: print "	Using ParMETIS"
else: print "	Not using ParMETIS (requires MPI)"
if env['usepapi']: print "	Using PAPI"
else: print "	Not using PAPI"
if env['usedebug']: print "	Compiling for debug"
else: print "	Not compiling for debug"
print "	Installing in", prefix
if ((fatalwarning != "") and (env['usewarnings'])): print "	Treating warnings as errors"
else: print "	Not treating warnings as errors"
print ""

############ Delete option-dependent files #####################

Execute(Delete(os.path.join(env['libinstall'],"Compiled.with.debug")))
Execute(Delete(os.path.join(env['libinstall'],"Compiled.with.mpi")))
Execute(Delete(os.path.join(env['libinstall'],"Compiled.with.openmp")))
Execute(Delete(os.path.join(env['libinstall'],"pyversion")))
Execute(Delete(os.path.join(env['libinstall'],"buildvars")))
if not env['usempi']: Execute(Delete(os.path.join(env['libinstall'],"pythonMPI")))


############ Build the subdirectories ##########################

from grouptest import *

TestGroups=[]

Export(
  ["env",
   "env_mpi",
   "clone_env",
   "IS_WINDOWS_PLATFORM",
   "TestGroups"
   ]
  )

env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
env.SConscript(dirs = ['tools/libescriptreader/src'], build_dir='build/$PLATFORM/tools/libescriptreader', duplicate=0)
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
env.SConscript(dirs = ['escript/src'], build_dir='build/$PLATFORM/escript', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['finley/src'], build_dir='build/$PLATFORM/finley', duplicate=0)
env.SConscript(dirs = ['modellib/py_src'], build_dir='build/$PLATFORM/modellib', duplicate=0)
env.SConscript(dirs = ['doc'], build_dir='build/$PLATFORM/doc', duplicate=0)
env.SConscript(dirs = ['pyvisi/py_src'], build_dir='build/$PLATFORM/pyvisi', duplicate=0)
env.SConscript(dirs = ['pycad/py_src'], build_dir='build/$PLATFORM/pycad', duplicate=0)
env.SConscript(dirs = ['pythonMPI/src'], build_dir='build/$PLATFORM/pythonMPI', duplicate=0)
env.SConscript(dirs = ['scripts'], build_dir='build/$PLATFORM/scripts', duplicate=0)
env.SConscript(dirs = ['paso/profiling'], build_dir='build/$PLATFORM/paso/profiling', duplicate=0)


############ Remember what optimizations we used ###############

remember_list = []

if env['usedebug']:
  remember_list += env.Command(os.path.join(env['libinstall'],"Compiled.with.debug"), None, Touch('$TARGET'))

if env['usempi']:
  remember_list += env.Command(os.path.join(env['libinstall'],"Compiled.with.mpi"), None, Touch('$TARGET'))

if env['useopenmp']:
  remember_list += env.Command(os.path.join(env['libinstall'],"Compiled.with.openmp"), None, Touch('$TARGET'))

env.Alias('remember_options', remember_list)


############### Record python interpreter version ##############

if not IS_WINDOWS_PLATFORM:
  versionstring="Python "+str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2])
  os.system("echo "+versionstring+" > "+os.path.join(env['libinstall'],"pyversion"))

############## Populate the buildvars file #####################

buildvars=open(os.path.join(env['libinstall'],'buildvars'),'w')
buildvars.write('python='+str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2])+'\n')

# Find the boost version by extracting it from version.hpp
boosthpp=open(os.path.join(env['boost_path'],'boost','version.hpp'))
boostversion='unknown'
try:
    for line in boosthpp:
        ver=re.match(r'#define BOOST_VERSION (\d+)',line)
        if ver:
            boostversion=ver.group(1)
except StopIteration: 
    pass
buildvars.write("boost="+boostversion+"\n")
buildvars.write("svn_revision="+str(global_revision)+"\n")
out="usedebug="
if env['usedebug']:
    out+="y"
else:
    out+="n"
out+="\nusempi="
if env['usempi']:
    out+="y"
else:
    out+="n"
out+="\nuseopenmp="
if env['useopenmp']:
    out+="y"
else:
    out+="n"
buildvars.write(out+"\n")
buildvars.write("mpi_flavour="+env['mpi_flavour']+'\n')

buildvars.close()


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
if env['usempi']:		build_all_list += ['target_pythonMPI_exe']
#if not IS_WINDOWS_PLATFORM:	build_all_list += ['target_escript_wrapper']
if env['usesilo']:	build_all_list += ['target_escript2silo']
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
if env['usempi']:		install_all_list += ['target_install_pythonMPI_exe']
#if not IS_WINDOWS_PLATFORM:	install_all_list += ['target_install_escript_wrapper']
if env['usesilo']:	install_all_list += ['target_install_escript2silo']
install_all_list += ['remember_options']
env.Alias('install_all', install_all_list)

# Default target is install
env.Default('install_all')

############ Targets to build and run the test suite ###########

env.Alias('build_cppunittest', ['target_install_cppunittest_headers', 'target_cppunittest_a'])
env.Alias('install_cppunittest', ['build_cppunittest', 'target_install_cppunittest_a'])
env.Alias('run_tests', ['install_all', 'target_install_cppunittest_a'])
env.Alias('all_tests', ['install_all', 'target_install_cppunittest_a', 'run_tests', 'py_tests'])
env.Alias('build_full',['install_all','build_tests','build_py_tests'])

############ Targets to build the documentation ################

env.Alias('docs', ['examples_tarfile', 'examples_zipfile', 'api_epydoc', 'api_doxygen', 'guide_pdf', 'guide_html','install_pdf'])

if not IS_WINDOWS_PLATFORM:
   try:
   	utest=open("utest.sh","w")
	build_platform=os.name		#Sometimes Mac python says it is posix
	if (build_platform=='posix') and platform.system()=="Darwin":
		build_platform='darwin'
	utest.write(GroupTest.makeHeader(build_platform))
	for tests in TestGroups:
	    utest.write(tests.makeString())
	utest.close()
	os.chmod("utest.sh",stat.S_IRWXU|stat.S_IRGRP|stat.S_IXGRP|stat.S_IROTH|stat.S_IXOTH)
	print "utest.sh written"
   except IOError:
	print "Error attempting to write unittests file."
	sys.exit(1)

