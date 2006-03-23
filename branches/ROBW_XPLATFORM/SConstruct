# top-level Scons configuration file for all esys13 modules
# Begin initialisation Section
# all of this section just intialises default environments and helper 
# scripts. You shouldn't need to modify this section.
EnsureSConsVersion(0,96,91)
EnsurePythonVersion(2,3)

import sys, os
# Add our extensions
if sys.path.count('scons')==0: sys.path.append('scons')
import scons_extensions

# Default options and options help text
# These are defaults and can be overridden using command line arguments or an options file.
# if the options_file or ARGUMENTS do not exist then the ones listed as default here are used
# DO NOT CHANGE THEM HERE
if ARGUMENTS.get('options_file',0):
   options_file = ARGUMENTS.get('options_file',0)
else:
   import socket
   from string import ascii_letters,digits
   hostname=""
   for s in socket.gethostname().split('.')[0]:
      if s in ascii_letters+digits:
         hostname+=s
      else:
         hostname+="_"
   options_file = "scons/"+hostname+"_options.py"

opts = Options(options_file, ARGUMENTS)
opts.AddOptions(
# Where to install esys stuff
  ('incinstall', 'where the esys headers will be installed', Dir('#.').abspath+'/include'), 
  ('libinstall', 'where the esys libraries will be installed', Dir('#.').abspath+'/lib'), 
  ('pyinstall', 'where the esys python modules will be installed', Dir('#.').abspath+'/esys'), 
# Compilation options
  BoolOption('dodebug', 'Do you want a debug build?', 'yes'),
  ('options_file', "Optional file containing preferred options. Ignored if it doesn't exist (default: scons/hostname_options.py)", options_file),
  ('cc_defines','C/C++ defines to use', None),
  ('cc_flags','C compiler flags to use (Release build)', None),
  ('cc_flags_debug', 'C compiler flags to use (Debug build)', None),
  ('cxx_flags', 'C++ compiler flags to use (Release build)', None),
  ('cxx_flags_debug', 'C++ compiler flags to use (Debug build)', None),
  ('ar_flags', 'Static library archiver flags to use', None),
  ('sys_libs', 'System libraries to link with', None),
# MKL
  PathOption('mkl_path', 'Path to MKL includes', None), 
  PathOption('mkl_lib_path', 'Path to MKL libs', None), 
  ('mkl_libs', 'MKL libraries to link with', None),
# SCSL
  PathOption('scsl_path', 'Path to SCSL includes', None), 
  PathOption('scsl_lib_path', 'Path to SCSL libs', None), 
  ('scsl_libs', 'SCSL libraries to link with', None),
# UMFPACK 
  PathOption('umf_path', 'Path to UMF includes', None), 
  PathOption('umf_lib_path', 'Path to UMF libs', None), 
  ('umf_libs', 'UMF libraries to link with', None),
# Python
  PathOption('python_path', 'Path to Python includes', '/usr/include'), 
  PathOption('python_lib_path', 'Path to Python libs', '/usr/lib'), 
  ('python_lib', 'Python libraries to link with', ['python2.3',]),
# Boost
  PathOption('boost_path', 'Path to Boost includes', '/usr/include'), 
  PathOption('boost_lib_path', 'Path to Boost libs', '/usr/lib'), 
  ('boost_lib', 'Boost libraries to link with', ['boost_python',]),
# CppUnit
  PathOption('cppunit_path', 'Path to CppUnit includes', Dir('#.').abspath+'/tools/CppUnitTest/inc'), 
  PathOption('cppunit_lib_path', 'Path to CppUnit libs', Dir('#.').abspath+'/lib'), 
  ('cppunit_lib', 'CppUnit libraries to link with', ['CppUnitTest',]),
# Doc building
  PathOption('doxygen_path', 'Path to Doxygen executable', None), 
  PathOption('epydoc_path', 'Path to Epydoc executable', None), 
  PathOption('epydoc_pythonpath', 'Path to Epydoc python files', None), 
# PAPI
  PathOption('papi_path', 'Path to PAPI includes', None), 
  PathOption('papi_lib_path', 'Path to PAPI libs', None), 
  ('papi_libs', 'PAPI libraries to link with', None),
)

# Initialise Scons Build Environment
# Note: On the Altix the intel compilers are not automatically
# detected by scons intelc.py script. The Altix has a different directory
# path and in some locations the "modules" facility is used to support
# multiple compiler versions. This forces the need to import the users PATH
# environment which isn't the "scons way"
# This doesn't impact linux and windows which will use the default compiler (g++ or msvc, or the intel compiler if it is installed on both platforms)
# FIXME: Perhaps a modification to intelc.py will allow better support for ia64 on altix
if os.name != "nt" and os.uname()[4]=='ia64':
   env = Environment(ENV = {'PATH' : os.environ['PATH']}, tools = ['default', 'intelc'], options = opts)
   if env['CXX'] == 'icpc':
      env['LINK'] = env['CXX'] # version >=9 of intel c++ compiler requires use of icpc to link in C++ runtimes (icc does not). FIXME: this behaviour could be directly incorporated into scons intelc.py
elif os.name == "nt":
   env = Environment(tools = ['default', 'intelc'], options = opts)
else:
   env = Environment(tools = ['default'], options = opts)

# Setup help for options
Help(opts.GenerateHelpText(env))

# Add some customer builders
py_builder = Builder(action = scons_extensions.build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

if env['PLATFORM'] == "win32":
   runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed', src_suffix='.exe', single_source=True)
else:
   runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed', single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

# FIXME: Copy this builder from the win32 branch
#runPyUnitTest_builder = Builder(action = scons_extensions.runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
#env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

# Convert the options which are held in environment variable into python variables for ease of handling and configure compilation options
try:
   incinstall = env['incinstall']
   env.Append(CPPPATH = [incinstall,])
except KeyError:
   incinstall = None  
try:
   libinstall = env['libinstall']
   env.Append(LIBPATH = [libinstall,])
except KeyError:
   libinstall = None  
try:
   pyinstall = env['pyinstall']
except KeyError:
   pyinstall = None  
try:
   dodebug = env['dodebug']
except KeyError:
   dodebug = None  
try:
   cc_defines = env['cc_defines']
   env.Append(CPPDEFINES = [cc_defines,])
except KeyError:
   pass
if dodebug:
   try:
      flags = env['cc_flags_debug']
      env.Append(CCFLAGS = flags)
   except KeyError:
      pass
else:
   try:
      flags = env['cc_flags']
      env.Append(CCFLAGS = flags)
   except KeyError:
      pass

if dodebug:
   try:
      flags = env['cxx_flags_debug']
      env.Append(CXXFLAGS = flags)
   except KeyError:
      pass
else:
   try:
      flags = env['cxx_flags']
      env.Append(CXXFLAGS = flags)
   except KeyError:
      pass

try:
   flags = env['ar_flags']
   env.Append(ARFLAGS = flags)
except KeyError:
   ar_flags = None  
try:
   sys_libs = env['sys_libs']
except KeyError:
   sys_libs = '' 

try:
   includes = env['mkl_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['mkl_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   mkl_libs = env['mkl_libs']
except KeyError:
   mkl_libs = ''
try:
   includes = env['scsl_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['scsl_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   scsl_libs = env['scsl_libs']
except KeyError:
   scsl_libs = ''
try:
   includes = env['umf_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['umf_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   umf_libs = env['umf_libs']
except KeyError:
   umf_libs = ''
try:
   includes = env['boost_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['boost_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   boost_lib = env['boost_lib']
except KeyError:
   boost_lib = None  
try:
   includes = env['cppunit_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['cppunit_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   cppunit_lib = env['cppunit_lib']
except KeyError:
   boost_lib = None  
try:
   includes = env['python_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['python_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   python_lib = env['python_lib']
except KeyError:
   python_lib = None  
try:
   doxygen_path = env['doxygen_path']
except KeyError:
   doxygen_path = None  
try:
   epydoc_path = env['epydoc_path']
except KeyError:
   epydoc_path = None  
try:
   epydoc_pythonpath = env['epydoc_pythonpath']
except KeyError:
   epydoc_pythonpath = None  
try:
   includes = env['papi_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['papi_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   papi_libs = env['papi_libs']
except KeyError:
   papi_libs = None  

# Targets
env.Default(libinstall)
env.Default(incinstall)
env.Alias('py_test')

# Allow sconscripts to see the env
Export(["env", "incinstall", "libinstall", "pyinstall", "dodebug", "mkl_libs", "scsl_libs", "umf_libs",
	"boost_lib", "python_lib", "doxygen_path", "epydoc_path", "epydoc_pythonpath", "papi_libs", "cppunit_lib", "sys_libs" ])

# End initialisation section
# Begin configuration section
# Insert new components to be build here
# FIXME: might be nice to replace this verbosity with a list of targets and some
# FIXME: nifty python to create the lengthy but very similar env.Sconscript lines
# Third Party libraries
env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
# C/C++ Libraries
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
env.SConscript(dirs = ['bruce/src'], build_dir='build/$PLATFORM/bruce', duplicate=0)
env.SConscript(dirs = ['escript/src'], build_dir='build/$PLATFORM/escript', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['finley/src'], build_dir='build/$PLATFORM/finley', duplicate=0)

# Unit Tests
#env.SConscript(dirs = ['bruce/test'], build_dir='build/$PLATFORM/bruce/test', duplicate=0)
#env.SConscript(dirs = ['finley/test'], build_dir='build/$PLATFORM/finley/test', duplicate=0)

# Python
#env.SConscript(dirs = ['esys/py_src'], build_dir='build/$PLATFORM/esys/py', duplicate=0)# call appropriate SConscripts
# FIXME:modelib and pyvisi need to be incorporated into build system for it to match original one
#'modellib/SConstruct',
#'pyvisi/SConstruct']
# 'doc/SConstruct']
# 'doc/SConstruct']

# FIXME: REMOVE THIS OLD STUFF
# Original SConstruct file contents (minus the bits that have been replaced)
#target_scripts = [
#                  'modellib/SConstruct',
#                  'pyvisi/SConstruct']
                  # 'doc/SConstruct']
                  # 'doc/SConstruct']

#SConscript(target_scripts, duplicate=0)

