opts = Options('custom.py')
opts.AddOptions(
   BoolOption('RELEASE', 'Set to build for release', 0),
   PathOption('PYTHON_HOME','Path to python home','C:/python23')
)
env = Environment(tools = ['default'],options = opts)

print "PLATFORM is:", env['PLATFORM']

EnsurePythonVersion(2,3)

#TODO: How do I convert these to options?
#TODO: Is there a more compact way of setting up the include paths? 

# Third-Party libraries 
boost_home = 'E:/woo409/development/boost'
python_home = env['PYTHON_HOME'] 

# Where to install (and find) esys includes and libraries
# Note: #/ means relative to the top of source tree
esys_inc = '#/include'
esys_lib = '#/lib'

# Derived paths
python_inc = python_home + '/include'
python_lib = python_home + '/libs'
boost_inc = boost_home
boost_lib = boost_home + '/windows_binary/lib'

incdir = [ boost_inc, python_inc, esys_inc ]
libdir = [ boost_lib, python_lib, esys_lib ]

env.Append(CPPPATH=incdir)
env.Append(LIBPATH=libdir)

env.Append(CPPDEFINES = ['DOASSERT'])

if env['PLATFORM'] == "win32":
   env.Append(CCFLAGS = ' /GR /EHsc /TP /wd4068')
   env.Append(CPPDEFINES = ['MSVC', 'WIN32'])
   if False :
      print "RELEASE build"
   else:
      print "DEBUG build"
      env.Append(CCFLAGS = ' /Od /MDd /RTC1')
      env.Append(CPPDEFINES = ['_DEBUG'])
      boost_lib_name = 'boost_python-vc71-mt-sgd'

Export(["env", "incdir", "esys_inc", "esys_lib", "boost_lib_name" ])

# Libraries
env.SConscript(dirs = ['paso/src'], build_dir='build/win32/paso', duplicate=0)
env.SConscript(dirs = ['bruce/src'], build_dir='build/win32/bruce', duplicate=0)
env.SConscript(dirs = ['escript/src/Data'], build_dir='build/win32/escript/Data', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/win32/esysUtils', duplicate=0)
env.SConscript(dirs = ['win32/win32_utils'], build_dir='build/win32/win32_utils', duplicate=0)
env.SConscript(dirs = ['tools/mmio/src'], build_dir='build/win32/tools/mmio', duplicate=0)
env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/win32/tools/CppUnitTest', duplicate=0)
env.SConscript(dirs = ['finley/src/finleyC'], build_dir='build/win32/finleyC', duplicate=0)
env.SConscript(dirs = ['finley/src/CPPAdapter'], build_dir='build/win32/CPPAdapter', duplicate=0)

# Unit Tests
env.SConscript(dirs = ['esysUtils/test/EsysException'], build_dir='build/win32/esysUtils/test/EsysException', duplicate=0)
env.SConscript(dirs = ['escript/test'], build_dir='build/win32/escript/test', duplicate=0)
env.SConscript(dirs = ['bruce/test'], build_dir='build/win32/bruce/test', duplicate=0)
env.SConscript(dirs = ['finley/test'], build_dir='build/win32/finley/test', duplicate=0)
