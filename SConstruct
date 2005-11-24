opts = Options('custom.py')
opts.AddOptions(
   BoolOption('RELEASE', 'Set to build for release', 0),
   PathOption('PYTHON_HOME','Path to python home','/usr/lib/python2.3')
)


# Extensions to Scons
def build_py(target, source, env):
   # Code to build .pyc from .py
   import py_compile, sys;
   py_compile.compile(str(source[0]), str(target[0]))
   return None
py_builder = Builder(action = build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)

env = Environment(tools = ['default'],options = opts)
env.Append(BUILDERS = {'PyCompile' : py_builder});

print "PLATFORM is:", env['PLATFORM']

EnsurePythonVersion(2,3)

#TODO: How do I convert these to options?
#TODO: Is there a more compact way of setting up the include paths? 

# Third-Party libraries 
boost_home = '/usr/include/boost'
python_home = env['PYTHON_HOME'] 

# Where to install (and find) esys includes and libraries
# Note: #/ means relative to the top of source tree
esys_inc = '#/include'
esys_lib = '#/lib'

# Derived paths
if env['PLATFORM'] == "win32":
   python_inc = python_home + '/include'
   python_lib = python_home + '/libs'
   boost_inc = boost_home
   boost_lib = boost_home + '/usr/lib'
elif  env['PLATFORM'] == "posix":
   python_inc = '/usr/include/python2.3'
   python_lib = '/usr/lib'
   boost_inc = '/usr/include'
   boost_lib = '/usr/lib'

incdir = [ boost_inc, python_inc, esys_inc ]
libdir = [ boost_lib, python_lib, esys_lib ]

env.Append(CPPPATH=incdir)
env.Append(LIBPATH=libdir)


if env['PLATFORM'] == "win32":
   env.Append(CCFLAGS = ' /GR /EHsc /TP /wd4068')
   env.Append(CPPDEFINES = ['MSVC', 'WIN32'])
   if False :
      print "RELEASE build"
   else:
      print "DEBUG build"
      env.Append(CPPDEFINES = ['DOASSERT'])
      env.Append(CCFLAGS = ' /Od /MDd /RTC1')
      env.Append(CPPDEFINES = ['_DEBUG'])
      boost_lib_name = 'boost_python-vc71-mt-sgd'
elif env['PLATFORM'] == "posix":
   env.Append(CC = ' -std=c99')
   env.Append(CCFLAGS = ' -c -fpic -W -Wall -Wno-unknown-pragmas')
   boost_lib_name = 'boost_python'
   if False :
      print "RELEASE build"
   else:
      print "DEBUG build"
      env.Append(CPPDEFINES = ['DOASSERT', 'DOPROF'])
      env.Prepend(CCFLAGS = ' -g -O0')

Export(["env", "incdir", "esys_inc", "esys_lib", "boost_lib_name" ])

# C/C++ Libraries
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
env.SConscript(dirs = ['bruce/src'], build_dir='build/$PLATFORM/bruce', duplicate=0)
env.SConscript(dirs = ['escript/src/Data'], build_dir='build/$PLATFORM/escript/Data', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['tools/mmio/src'], build_dir='build/$PLATFORM/tools/mmio', duplicate=0)
env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
env.SConscript(dirs = ['finley/src/finleyC'], build_dir='build/$PLATFORM/finleyC', duplicate=0)
env.SConscript(dirs = ['finley/src/CPPAdapter'], build_dir='build/$PLATFORM/CPPAdapter', duplicate=0)

# Unit Tests
env.SConscript(dirs = ['esysUtils/test/EsysException'], build_dir='build/$PLATFORM/esysUtils/test/EsysException', duplicate=0)
env.SConscript(dirs = ['escript/test'], build_dir='build/$PLATFORM/escript/test', duplicate=0)
env.SConscript(dirs = ['bruce/test'], build_dir='build/$PLATFORM/bruce/test', duplicate=0)
env.SConscript(dirs = ['finley/test'], build_dir='build/$PLATFORM/finley/test', duplicate=0)

if env['PLATFORM'] == "win32":
   env.SConscript(dirs = ['win32/win32_utils'], build_dir='build/$PLATFORM/win32_utils', duplicate=0)
