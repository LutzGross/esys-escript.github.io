import scons_ext

opts = Options('custom.py')
opts.AddOptions(
    BoolOption('RELEASE', 'Set to build for release', False),
    PathOption('PYTHON_INC','Path to python includes','/usr/include/python2.3'),
    PathOption('PYTHON_LIB','Path to python includes','/usr/lib'),
    PathOption('BOOST_INC','Path to boost includes','/usr/include'),
    PathOption('BOOST_LIB','Path to boost includes','/usr/lib'),

    # CC and CXX overrides of default environments
    ('CC_PATH','override CC',False),
    ('CXX_PATH','override CXX',False)
)

env = Environment(tools = ['default'],options = opts)

if env['CC_PATH'] :
    env['CC'] = env['CC_PATH']

if env['CXX_PATH'] :
    env['CXX'] = env['CXX_PATH']

Help(opts.GenerateHelpText(env))

py_builder = Builder(action = scons_ext.build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

if env['PLATFORM'] == "win32":
    runUnitTest_builder = Builder(action = scons_ext.runUnitTest, suffix = '.passed', src_suffix='.exe', single_source=True)
else:
    runUnitTest_builder = Builder(action = scons_ext.runUnitTest, suffix = '.passed', single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = scons_ext.runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});


print "PLATFORM is:", env['PLATFORM']

EnsurePythonVersion(2,3)

# Third-Party libraries 
python_inc = env['PYTHON_INC']
python_lib = env['PYTHON_LIB']
boost_inc = env['BOOST_INC']
boost_lib = env['BOOST_LIB']

# Where to install (and find) esys includes and libraries
# Note: #/ means relative to the top of source tree
esys_inc = '#/include'
esys_lib = '#/lib'

env.Default(esys_lib)
env.Alias('py_test')

incdir = [ boost_inc, python_inc, esys_inc ]
libdir = [ boost_lib, python_lib, esys_lib ]

env.Append(CPPPATH=incdir)
env.Append(LIBPATH=libdir)

if not env['RELEASE'] :
    env.Append(CPPDEFINES = ['DOASSERT' 'DOPROF'])

if  env['CC'].find('cl.exe') >= 0 :
    env.Append(CCFLAGS = ' /GR /EHsc /TP /wd4068')
    env.Append(CPPDEFINES = ['MSVC', 'WIN32'])
    if env['RELEASE'] :
        print "RELEASE build"
    else:
       print "DEBUG build"
       env.Append(CCFLAGS = ' /Od /MDd /RTC1')
       env.Append(CPPDEFINES = ['_DEBUG'])
       boost_lib_name = 'boost_python-vc71-mt-sgd'

elif env['CC'].find('gcc') >= 0 :
    env.Append(CC = ' -std=c99')
    env.Append(CCFLAGS = ' -c -fpic -W -Wall -Wno-unknown-pragmas')
    boost_lib_name = 'boost_python'
    if env['RELEASE'] :
        print "RELEASE build"
        env.Prepend(CCFLAGS = ' -O3')
    else:
        print "DEBUG build"
        env.Prepend(CCFLAGS = ' -g -O0')

elif env['CC'].find('icc') >= 0 :
    env.Append(CC = ' -c99')
    env.Append(CXX = ' -ansi')
    env.Append(CCFLAGS = '  -openmp -openmp_report0 -mp1 -tpp2 -ansi_alias -no-gcc -fpic -w1')
    env.Append(CPPDEFINES = ['SCSL'])
    boost_lib_name = 'boost_python'
    if env['RELEASE'] :
        print "RELEASE build"
        env.Prepend(CCFLAGS = ' -O3 -IPF_fma -ftz')
    else:
        print "DEBUG build"
        env.Prepend(CCFLAGS = ' -g -O0')

Export(["env", "incdir", "esys_inc", "esys_lib", "boost_lib_name" ])

# C/C++ Libraries
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
env.SConscript(dirs = ['bruce/src'], build_dir='build/$PLATFORM/bruce', duplicate=0)
env.SConscript(dirs = ['escript/src/Data'], build_dir='build/$PLATFORM/escript/Data', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['tools/mmio/src'], build_dir='build/$PLATFORM/tools/mmio', duplicate=0)
env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
env.SConscript(dirs = ['finley/src/finleyC'], build_dir='build/$PLATFORM/finley/finleyC', duplicate=0)
env.SConscript(dirs = ['finley/src/CPPAdapter'], build_dir='build/$PLATFORM/finley/CPPAdapter', duplicate=0)

if env['PLATFORM'] == "win32":
    env.SConscript(dirs = ['win32/win32_utils'], build_dir='build/$PLATFORM/win32_utils', duplicate=0)

# Unit Tests
env.SConscript(dirs = ['esysUtils/test/EsysException'], build_dir='build/$PLATFORM/esysUtils/test/EsysException', duplicate=0)
env.SConscript(dirs = ['escript/test'], build_dir='build/$PLATFORM/escript/test', duplicate=0)
env.SConscript(dirs = ['bruce/test'], build_dir='build/$PLATFORM/bruce/test', duplicate=0)
env.SConscript(dirs = ['finley/test'], build_dir='build/$PLATFORM/finley/test', duplicate=0)

# Python
env.SConscript(dirs = ['esys/py_src'], build_dir='build/$PLATFORM/esys/py', duplicate=0)
