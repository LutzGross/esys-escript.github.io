##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

EnsureSConsVersion(0,98,1)
EnsurePythonVersion(2,5)

import atexit, sys, os, platform, re
from distutils import sysconfig
from dependencies import *
from site_init import *

# Version number to check for in options file. Increment when new features are
# added or existing options changed.
REQUIRED_OPTS_VERSION=201

# MS Windows support, many thanks to PH
IS_WINDOWS = (os.name == 'nt')

########################## Determine options file ############################
# 1. command line
# 2. scons/<hostname>_options.py
# 3. name as part of a cluster
options_file=ARGUMENTS.get('options_file', None)
if not options_file:
    ext_dir = os.path.join(os.getcwd(), 'scons')
    hostname = platform.node().split('.')[0]
    for name in hostname, effectiveName(hostname):
        mangledhostname = re.sub('[^0-9a-zA-Z]', '_', hostname)
        options_file = os.path.join(ext_dir, mangledhostname+'_options.py')
        if os.path.isfile(options_file): break

if not os.path.isfile(options_file):
    print("\nWARNING:\nOptions file %s" % options_file)
    print("not found! Default options will be used which is most likely suboptimal.")
    print("It is recommended that you copy one of the TEMPLATE files in the scons/")
    print("subdirectory and customize it to your needs.\n")
    options_file = None

############################### Build options ################################

default_prefix='/usr'
mpi_flavours=('no', 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI')
lapack_flavours=('none', 'clapack', 'mkl')

vars = Variables(options_file, ARGUMENTS)
vars.AddVariables(
  PathVariable('options_file', 'Path to options file', options_file, PathVariable.PathIsFile),
  PathVariable('prefix', 'Installation prefix', Dir('#.').abspath, PathVariable.PathIsDirCreate),
  PathVariable('build_dir', 'Top-level build directory', Dir('#/build').abspath, PathVariable.PathIsDirCreate),
  BoolVariable('verbose', 'Output full compile/link lines', False),
# Compiler/Linker options
  ('cc', 'Path to C compiler', 'default'),
  ('cxx', 'Path to C++ compiler', 'default'),
  ('cc_flags', 'Base C/C++ compiler flags', 'default'),
  ('cc_optim', 'Additional C/C++ flags for a non-debug build', 'default'),
  ('cc_debug', 'Additional C/C++ flags for a debug build', 'default'),
  ('cc_extra', 'Extra C compiler flags', ''),
  ('cxx_extra', 'Extra C++ compiler flags', ''),
  ('ld_extra', 'Extra linker flags', ''),
  BoolVariable('werror','Treat compiler warnings as errors', True),
  BoolVariable('debug', 'Compile with debug flags', False),
  BoolVariable('openmp', 'Compile parallel version using OpenMP', False),
  ('omp_flags', 'OpenMP compiler flags', 'default'),
  ('omp_ldflags', 'OpenMP linker flags', 'default'),
# Mandatory libraries
  ('boost_prefix', 'Prefix/Paths of boost installation', default_prefix),
  ('boost_libs', 'Boost libraries to link with', ['boost_python-mt']),
# Mandatory for tests
  ('cppunit_prefix', 'Prefix/Paths of CppUnit installation', default_prefix),
  ('cppunit_libs', 'CppUnit libraries to link with', ['cppunit']),
# Optional libraries and options
  EnumVariable('mpi', 'Compile parallel version using MPI flavour', 'none', allowed_values=mpi_flavours),
  ('mpi_prefix', 'Prefix/Paths of MPI installation', default_prefix),
  ('mpi_libs', 'MPI shared libraries to link with', ['mpi']),
  BoolVariable('netcdf', 'Enable netCDF file support', False),
  ('netcdf_prefix', 'Prefix/Paths of netCDF installation', default_prefix),
  ('netcdf_libs', 'netCDF libraries to link with', ['netcdf_c++', 'netcdf']),
  BoolVariable('parmetis', 'Enable ParMETIS (requires MPI)', False),
  ('parmetis_prefix', 'Prefix/Paths of ParMETIS installation', default_prefix),
  ('parmetis_libs', 'ParMETIS libraries to link with', ['parmetis', 'metis']),
  BoolVariable('papi', 'Enable PAPI', False),
  ('papi_prefix', 'Prefix/Paths to PAPI installation', default_prefix),
  ('papi_libs', 'PAPI libraries to link with', ['papi']),
  BoolVariable('papi_instrument_solver', 'Use PAPI to instrument each iteration of the solver', False),
  BoolVariable('mkl', 'Enable the Math Kernel Library', False),
  ('mkl_prefix', 'Prefix/Paths to MKL installation', default_prefix),
  ('mkl_libs', 'MKL libraries to link with', ['mkl_solver','mkl_em64t','guide','pthread']),
  BoolVariable('umfpack', 'Enable UMFPACK', False),
  ('umfpack_prefix', 'Prefix/Paths to UMFPACK installation', default_prefix),
  ('umfpack_libs', 'UMFPACK libraries to link with', ['umfpack']),
  BoolVariable('boomeramg', 'Enable BoomerAMG', False),
  ('boomeramg_prefix', 'Prefix/Paths to BoomerAMG installation', default_prefix),
  ('boomeramg_libs', 'BoomerAMG libraries to link with', ['boomeramg']),
  EnumVariable('lapack', 'Set LAPACK flavour', 'none', allowed_values=lapack_flavours),
  ('lapack_prefix', 'Prefix/Paths to LAPACK installation', default_prefix),
  ('lapack_libs', 'LAPACK libraries to link with', []),
  BoolVariable('silo', 'Enable the Silo file format in weipa', False),
  ('silo_prefix', 'Prefix/Paths to Silo installation', default_prefix),
  ('silo_libs', 'Silo libraries to link with', ['siloh5', 'hdf5']),
  BoolVariable('visit', 'Enable the VisIt simulation interface', False),
  ('visit_prefix', 'Prefix/Paths to VisIt installation', default_prefix),
  ('visit_libs', 'VisIt libraries to link with', ['simV2']),
  BoolVariable('vsl_random', 'Use VSL from intel for random data', False),
# Advanced settings
  #dudley_assemble_flags = -funroll-loops      to actually do something
  ('dudley_assemble_flags', 'compiler flags for some dudley optimisations', ''),
  # To enable passing function pointers through python
  BoolVariable('iknowwhatimdoing', 'Allow non-standard C', False),
  # An option for specifying the compiler tools (see windows branch)
  ('tools_names', 'Compiler tools to use', ['default']),
  ('env_export', 'Environment variables to be passed to tools',[]),
  EnumVariable('forcelazy', 'For testing use only - set the default value for autolazy', 'leave_alone', allowed_values=('leave_alone', 'on', 'off')),
  EnumVariable('forcecollres', 'For testing use only - set the default value for force resolving collective ops', 'leave_alone', allowed_values=('leave_alone', 'on', 'off')),
  # finer control over library building, intel aggressive global optimisation
  # works with dynamic libraries on windows.
  ('build_shared', 'Build dynamic libraries only', False),
  ('sys_libs', 'Extra libraries to link with', []),
  ('escript_opts_version', 'Version of options file (do not specify on command line)'),
  ('SVN_VERSION', 'Do not use from options file', -2),
  ('pythoncmd', 'which python to compile with','python'),
  ('usepython3', 'Is this a python3 build? (experimental)', False),
  ('pythonlibname', 'Name of the python library to link. (This is found automatically for python2.X.)', ''),
  ('pythonlibpath', 'Path to the python library. (You should not need to set this unless your python has moved)',''),
  ('pythonincpath','Path to python include files. (You should not need to set this unless your python has moved',''),
  BoolVariable('BADPYTHONMACROS','Extra \#include to get around a python bug.', True),
)

##################### Create environment and help text #######################

# Intel's compiler uses regular expressions improperly and emits a warning
# about failing to find the compilers. This warning can be safely ignored.

# PATH is needed so the compiler, linker and tools are found if they are not
# in default locations.
env = Environment(tools = ['default'], options = vars,
                  ENV = {'PATH': os.environ['PATH']})

# set the vars for clang
def mkclang(env):
    env['CC']='clang'
    env['CXX']='clang++'

if env['tools_names'] != 'default':
    zz=env['tools_names']
    if 'clang' in zz:
        zz.remove('clang')
        zz.insert(0, mkclang)
    env = Environment(tools = ['default'] + env['tools_names'], options = vars,
                      ENV = {'PATH' : os.environ['PATH']})

if options_file:
    opts_valid=False
    if 'escript_opts_version' in env.Dictionary() and \
        int(env['escript_opts_version']) >= REQUIRED_OPTS_VERSION:
            opts_valid=True
    if opts_valid:
        print("Using options in %s." % options_file)
    else:
        print("\nOptions file %s" % options_file)
        print("is outdated! Please update the file by examining one of the TEMPLATE")
        print("files in the scons/ subdirectory and setting escript_opts_version to %d.\n"%REQUIRED_OPTS_VERSION)
        Exit(1)

# Generate help text (scons -h)
Help(vars.GenerateHelpText(env))

# Check for superfluous options
if len(vars.UnknownVariables())>0:
    for k in vars.UnknownVariables():
        print("Unknown option '%s'" % k)
    Exit(1)

# create dictionary which will be populated with info for buildvars file
env['buildvars']={}
# create list which will be populated with warnings if there are any
env['warnings']=[]

#################### Make sure install directories exist #####################

env['BUILD_DIR']=Dir(env['build_dir']).abspath
prefix=Dir(env['prefix']).abspath
env['buildvars']['prefix']=prefix
env['incinstall'] = os.path.join(prefix, 'include')
env['bininstall'] = os.path.join(prefix, 'bin')
env['libinstall'] = os.path.join(prefix, 'lib')
env['pyinstall']  = os.path.join(prefix, 'esys')
if not os.path.isdir(env['bininstall']):
    os.makedirs(env['bininstall'])
if not os.path.isdir(env['libinstall']):
    os.makedirs(env['libinstall'])
if not os.path.isdir(env['pyinstall']):
    os.makedirs(env['pyinstall'])

env.Append(CPPPATH = [env['incinstall']])
env.Append(LIBPATH = [env['libinstall']])

################# Fill in compiler options if not set above ##################

if env['cc'] != 'default': env['CC']=env['cc']
if env['cxx'] != 'default': env['CXX']=env['cxx']

# version >=9 of intel C++ compiler requires use of icpc to link in C++
# runtimes (icc does not)
if not IS_WINDOWS and os.uname()[4]=='ia64' and env['CXX']=='icpc':
    env['LINK'] = env['CXX']

# default compiler/linker options
cc_flags = ''
cc_optim = ''
cc_debug = ''
omp_flags = ''
omp_ldflags = ''
fatalwarning = '' # switch to turn warnings into errors
sysheaderopt = '' # how to indicate that a header is a system header

# env['CC'] might be a full path
cc_name=os.path.basename(env['CC'])

if cc_name == 'icc':
    # Intel compiler
    # #1875: offsetof applied to non-POD types is nonstandard (in boost)
    cc_flags    = "-std=c99 -fPIC -w2 -wd1875 -Wno-unknown-pragmas -DBLOCKTIMER -DCORE_ID1"
    cc_optim    = "-O3 -ftz -fno-alias -ipo -xHost"
    cc_debug    = "-g -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK"
    omp_flags   = "-openmp"
    omp_ldflags = "-openmp -openmp_report=1"
    fatalwarning = "-Werror"
elif cc_name[:3] == 'gcc':
    # GNU C on any system
    # note that -ffast-math is not used because it breaks isnan(),
    # see mantis #691
    cc_flags     = "-pedantic -Wall -fPIC -Wno-unknown-pragmas -DBLOCKTIMER  -Wno-sign-compare -Wno-system-headers -Wno-long-long -Wno-strict-aliasing -finline-functions"
    cc_optim     = "-O3"
    cc_debug     = "-g -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK"
    omp_flags    = "-fopenmp"
    omp_ldflags  = "-fopenmp"
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"
elif cc_name == 'cl':
    # Microsoft Visual C on Windows
    cc_flags     = "/EHsc /MD /GR /wd4068 /D_USE_MATH_DEFINES /DDLL_NETCDF"
    cc_optim     = "/O2 /Op /W3"
    cc_debug     = "/Od /RTCcsu /ZI /DBOUNDS_CHECK"
    fatalwarning = "/WX"
elif cc_name == 'icl':
    # Intel C on Windows
    cc_flags     = '/EHsc /GR /MD'
    cc_optim     = '/fast /Oi /W3 /Qssp /Qinline-factor- /Qinline-min-size=0 /Qunroll'
    cc_debug     = '/Od /RTCcsu /Zi /Y- /debug:all /Qtrapuv'
    omp_flags    = '/Qvec-report0 /Qopenmp /Qopenmp-report0 /Qparallel'
    omp_ldflags  = '/Qvec-report0 /Qopenmp /Qopenmp-report0 /Qparallel'

env['sysheaderopt']=sysheaderopt

# set defaults if not otherwise specified
if env['cc_flags']    == 'default': env['cc_flags'] = cc_flags
if env['cc_optim']    == 'default': env['cc_optim'] = cc_optim
if env['cc_debug']    == 'default': env['cc_debug'] = cc_debug
if env['omp_flags']   == 'default': env['omp_flags'] = omp_flags
if env['omp_ldflags'] == 'default': env['omp_ldflags'] = omp_ldflags
if env['cc_extra']  != '': env.Append(CFLAGS = env['cc_extra'])
if env['cxx_extra'] != '': env.Append(CXXFLAGS = env['cxx_extra'])
if env['ld_extra']  != '': env.Append(LINKFLAGS = env['ld_extra'])

if env['BADPYTHONMACROS']: env.Append(CXXFLAGS = ' -DBADPYTHONMACROS')

if env['usepython3']:
    env.Append(CPPDEFINES=['ESPYTHON3'])

# set up the autolazy values
if env['forcelazy'] == 'on':
    env.Append(CPPDEFINES=['FAUTOLAZYON'])
elif env['forcelazy'] == 'off':
    env.Append(CPPDEFINES=['FAUTOLAZYOFF'])

# set up the collective resolve values
if env['forcecollres'] == 'on':
    env.Append(CPPDEFINES=['FRESCOLLECTON'])
elif env['forcecollres'] == 'off':
    env.Append(CPPDEFINES=['FRESCOLLECTOFF'])

# allow non-standard C if requested
if env['iknowwhatimdoing']:
    env.Append(CPPDEFINES=['IKNOWWHATIMDOING'])

# Disable OpenMP if no flags provided
if env['openmp'] and env['omp_flags'] == '':
   env['warnings'].append("OpenMP requested but no flags provided - disabling OpenMP!")
   env['openmp'] = False

if env['openmp']:
    env.Append(CCFLAGS = env['omp_flags'])
    if env['omp_ldflags'] != '': env.Append(LINKFLAGS = env['omp_ldflags'])
else:
    env['omp_flags']=''
    env['omp_ldflags']=''

env['buildvars']['openmp']=int(env['openmp'])

# add debug/non-debug compiler flags
env['buildvars']['debug']=int(env['debug'])
if env['debug']:
    env.Append(CCFLAGS = env['cc_debug'])
else:
    env.Append(CCFLAGS = env['cc_optim'])

# always add cc_flags
env.Append(CCFLAGS = env['cc_flags'])

# add system libraries
env.AppendUnique(LIBS = env['sys_libs'])

# determine svn revision
global_revision=ARGUMENTS.get('SVN_VERSION', None)
if global_revision:
    global_revision = re.sub(':.*', '', global_revision)
    global_revision = re.sub('[^0-9]', '', global_revision)
    if global_revision == '': global_revision='-2'
else:
  # Get the global Subversion revision number for the getVersion() method
  try:
    global_revision = os.popen('svnversion -n .').read()
    global_revision = re.sub(':.*', '', global_revision)
    global_revision = re.sub('[^0-9]', '', global_revision)
    if global_revision == '': global_revision='-2'
  except:
    global_revision = '-1'
env['svn_revision']=global_revision
env['buildvars']['svn_revision']=global_revision
env.Append(CPPDEFINES=['SVN_VERSION='+global_revision])

if IS_WINDOWS:
    if not env['build_shared']:
        env.Append(CPPDEFINES = ['ESYSUTILS_STATIC_LIB'])
        env.Append(CPPDEFINES = ['PASO_STATIC_LIB'])

# VSL random numbers
env['buildvars']['vsl_random']=int(env['vsl_random'])
if env['vsl_random']:
    env.Append(CPPDEFINES = ['MKLRANDOM'])

env['IS_WINDOWS']=IS_WINDOWS

###################### Copy required environment vars ########################

# Windows doesn't use LD_LIBRARY_PATH but PATH instead
if IS_WINDOWS:
    LD_LIBRARY_PATH_KEY='PATH'
    env['ENV']['LD_LIBRARY_PATH']=''
else:
    LD_LIBRARY_PATH_KEY='LD_LIBRARY_PATH'

env['LD_LIBRARY_PATH_KEY']=LD_LIBRARY_PATH_KEY

# the following env variables are exported for the unit tests

for key in 'OMP_NUM_THREADS', 'ESCRIPT_NUM_PROCS', 'ESCRIPT_NUM_NODES':
    try:
        env['ENV'][key] = os.environ[key]
    except KeyError:
        env['ENV'][key] = 1

env_export=env['env_export']
env_export.extend(['ESCRIPT_NUM_THREADS','ESCRIPT_HOSTFILE','DISPLAY','XAUTHORITY','PATH','HOME','KMP_MONITOR_STACKSIZE','TMPDIR','TEMP','TMP'])

for key in set(env_export):
    try:
        env['ENV'][key] = os.environ[key]
    except KeyError:
        pass

try:
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, os.environ[LD_LIBRARY_PATH_KEY])
except KeyError:
    pass

# these shouldn't be needed
#for key in 'C_INCLUDE_PATH','CPLUS_INCLUDE_PATH','LIBRARY_PATH':
#    try:
#        env['ENV'][key] = os.environ[key]
#    except KeyError:
#        pass

try:
    env['ENV']['PYTHONPATH'] = os.environ['PYTHONPATH']
except KeyError:
    pass

######################## Add some custom builders ############################

if env['pythoncmd']=='python':
    py_builder = Builder(action = build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
else:
    py_builder = Builder(action = env['pythoncmd']+" scripts/py_comp.py $SOURCE $TARGET", suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

runUnitTest_builder = Builder(action = runUnitTest, suffix = '.passed', src_suffix=env['PROGSUFFIX'], single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

epstopdfbuilder = Builder(action = eps2pdf, suffix='.pdf', src_suffix='.eps', single_source=True)
env.Append(BUILDERS = {'EpsToPDF' : epstopdfbuilder});

############################ Dependency checks ###############################

######## Compiler
env=checkCompiler(env)

######## Python headers & library (required)
env=checkPython(env)

######## boost & boost-python (required)
env=checkBoost(env)

######## numpy (required) and numpy headers (optional)
env=checkNumpy(env)

######## CppUnit (required for tests)
env=checkCppUnit(env)

######## optional python modules (sympy, pyproj)
env=checkOptionalModules(env)

######## optional dependencies (netCDF, PAPI, MKL, UMFPACK, Lapack, Silo, ...)
env=checkOptionalLibraries(env)

######## PDFLaTeX (for documentation)
env=checkPDFLatex(env)

# keep some of our install paths first in the list for the unit tests
env.PrependENVPath(LD_LIBRARY_PATH_KEY, env['libinstall'])
env.PrependENVPath('PYTHONPATH', prefix)
env['ENV']['ESCRIPT_ROOT'] = prefix

if not env['verbose']:
    env['CCCOMSTR'] = "Compiling $TARGET"
    env['CXXCOMSTR'] = "Compiling $TARGET"
    env['SHCCCOMSTR'] = "Compiling $TARGET"
    env['SHCXXCOMSTR'] = "Compiling $TARGET"
    env['ARCOMSTR'] = "Linking $TARGET"
    env['LINKCOMSTR'] = "Linking $TARGET"
    env['SHLINKCOMSTR'] = "Linking $TARGET"
    env['PDFLATEXCOMSTR'] = "Building $TARGET from LaTeX input $SOURCES"
    env['BIBTEXCOMSTR'] = "Generating bibliography $TARGET"
    env['MAKEINDEXCOMSTR'] = "Generating index $TARGET"
    env['PDFLATEXCOMSTR'] = "Building $TARGET from LaTeX input $SOURCES"
    #Progress(['Checking -\r', 'Checking \\\r', 'Checking |\r', 'Checking /\r'], interval=17)

####################### Configure the subdirectories #########################

# remove obsolete files
if not env['usempi']:
    Execute(Delete(os.path.join(env['libinstall'], 'pythonMPI')))
    Execute(Delete(os.path.join(env['libinstall'], 'pythonMPIredirect')))

from grouptest import *
TestGroups=[]

# keep an environment without warnings-as-errors
dodgy_env=env.Clone()

# now add warnings-as-errors flags. This needs to be done after configuration
# because the scons test files have warnings in them
if ((fatalwarning != '') and (env['werror'])):
    env.Append(CCFLAGS = fatalwarning)

Export(
  ['env',
   'dodgy_env',
   'IS_WINDOWS',
   'TestGroups'
  ]
)

env.SConscript(dirs = ['tools/escriptconvert'], variant_dir='$BUILD_DIR/$PLATFORM/tools/escriptconvert', duplicate=0)
env.SConscript(dirs = ['paso/src'], variant_dir='$BUILD_DIR/$PLATFORM/paso', duplicate=0)
env.SConscript(dirs = ['weipa/src'], variant_dir='$BUILD_DIR/$PLATFORM/weipa', duplicate=0)
env.SConscript(dirs = ['escript/src'], variant_dir='$BUILD_DIR/$PLATFORM/escript', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], variant_dir='$BUILD_DIR/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['pasowrap/src'], variant_dir='$BUILD_DIR/$PLATFORM/pasowrap', duplicate=0)
env.SConscript(dirs = ['dudley/src'], variant_dir='$BUILD_DIR/$PLATFORM/dudley', duplicate=0)
env.SConscript(dirs = ['finley/src'], variant_dir='$BUILD_DIR/$PLATFORM/finley', duplicate=0)
env.SConscript(dirs = ['ripley/src'], variant_dir='$BUILD_DIR/$PLATFORM/ripley', duplicate=0)
env.SConscript(dirs = ['downunder/py_src'], variant_dir='$BUILD_DIR/$PLATFORM/downunder', duplicate=0)
env.SConscript(dirs = ['modellib/py_src'], variant_dir='$BUILD_DIR/$PLATFORM/modellib', duplicate=0)
env.SConscript(dirs = ['pycad/py_src'], variant_dir='$BUILD_DIR/$PLATFORM/pycad', duplicate=0)
env.SConscript(dirs = ['pythonMPI/src'], variant_dir='$BUILD_DIR/$PLATFORM/pythonMPI', duplicate=0)
env.SConscript(dirs = ['doc'], variant_dir='$BUILD_DIR/$PLATFORM/doc', duplicate=0)
env.SConscript(dirs = ['paso/profiling'], variant_dir='$BUILD_DIR/$PLATFORM/paso/profiling', duplicate=0)


######################## Populate the buildvars file #########################

write_buildvars(env)

################### Targets to build and install libraries ###################

target_init = env.Command(os.path.join(env['pyinstall'],'__init__.py'), None, Touch('$TARGET'))
env.Alias('target_init', [target_init])
# delete buildvars upon cleanup
env.Clean('target_init', os.path.join(env['libinstall'], 'buildvars'))

# The headers have to be installed prior to build in order to satisfy
# #include <paso/Common.h>
env.Alias('build_esysUtils', ['install_esysUtils_headers', 'build_esysUtils_lib'])
env.Alias('install_esysUtils', ['build_esysUtils', 'install_esysUtils_lib'])

env.Alias('build_paso', ['install_paso_headers', 'build_paso_lib'])
env.Alias('install_paso', ['build_paso', 'install_paso_lib'])

env.Alias('build_escript', ['install_escript_headers', 'build_escript_lib', 'build_escriptcpp_lib'])
env.Alias('install_escript', ['build_escript', 'install_escript_lib', 'install_escriptcpp_lib', 'install_escript_py'])

env.Alias('build_pasowrap', ['install_pasowrap_headers', 'build_pasowrap_lib', 'build_pasowrapcpp_lib'])
env.Alias('install_pasowrap', ['build_pasowrap', 'install_pasowrap_lib', 'install_pasowrapcpp_lib', 'install_pasowrap_py'])

env.Alias('build_dudley', ['install_dudley_headers', 'build_dudley_lib', 'build_dudleycpp_lib'])
env.Alias('install_dudley', ['build_dudley', 'install_dudley_lib', 'install_dudleycpp_lib', 'install_dudley_py'])

env.Alias('build_finley', ['install_finley_headers', 'build_finley_lib', 'build_finleycpp_lib'])
env.Alias('install_finley', ['build_finley', 'install_finley_lib', 'install_finleycpp_lib', 'install_finley_py'])

env.Alias('build_ripley', ['install_ripley_headers', 'build_ripley_lib', 'build_ripleycpp_lib'])
env.Alias('install_ripley', ['build_ripley', 'install_ripley_lib', 'install_ripleycpp_lib', 'install_ripley_py'])

env.Alias('build_weipa', ['install_weipa_headers', 'build_weipa_lib', 'build_weipacpp_lib'])
env.Alias('install_weipa', ['build_weipa', 'install_weipa_lib', 'install_weipacpp_lib', 'install_weipa_py'])

env.Alias('build_escriptreader', ['install_weipa_headers', 'build_escriptreader_lib'])
env.Alias('install_escriptreader', ['build_escriptreader', 'install_escriptreader_lib'])

# Now gather all the above into some easy targets: build_all and install_all
build_all_list = []
build_all_list += ['build_esysUtils']
build_all_list += ['build_paso']
build_all_list += ['build_escript']
build_all_list += ['build_pasowrap']
build_all_list += ['build_dudley']
build_all_list += ['build_finley']
build_all_list += ['build_ripley']
build_all_list += ['build_weipa']
if not IS_WINDOWS: build_all_list += ['build_escriptreader']
if env['usempi']:   build_all_list += ['build_pythonMPI']
build_all_list += ['build_escriptconvert']
env.Alias('build_all', build_all_list)

install_all_list = []
install_all_list += ['target_init']
install_all_list += ['install_esysUtils']
install_all_list += ['install_paso']
install_all_list += ['install_escript']
install_all_list += ['install_pasowrap']
install_all_list += ['install_dudley']
install_all_list += ['install_finley']
install_all_list += ['install_ripley']
install_all_list += ['install_weipa']
if not IS_WINDOWS: install_all_list += ['install_escriptreader']
install_all_list += ['install_downunder_py']
install_all_list += ['install_modellib_py']
install_all_list += ['install_pycad_py']
if env['usempi']:   install_all_list += ['install_pythonMPI']
install_all_list += ['install_escriptconvert']
env.Alias('install_all', install_all_list)

# Default target is install
env.Default('install_all')

################## Targets to build and run the test suite ###################

if not env['cppunit']:
    test_msg = env.Command('.dummy.', None, '@echo "Cannot run C/C++ unit tests, CppUnit not found!";exit 1')
    env.Alias('run_tests', test_msg)
    env.Alias('build_tests', '')
env.Alias('run_tests', ['install_all'])
env.Alias('all_tests', ['install_all', 'run_tests', 'py_tests'])
env.Alias('build_full',['install_all','build_tests','build_py_tests'])
env.Alias('build_PasoTests','$BUILD_DIR/$PLATFORM/paso/profiling/PasoTests')

##################### Targets to build the documentation #####################

env.Alias('pdfdocs',['user_pdf', 'install_pdf', 'cookbook_pdf', 'inversion_pdf'])
env.Alias('basedocs', ['pdfdocs','examples_tarfile', 'examples_zipfile', 'api_doxygen'])
env.Alias('docs', ['basedocs', 'sphinxdoc'])
env.Alias('release_prep', ['docs', 'install_all'])
env.Alias('release_prep_old', ['basedocs', 'api_epydoc', 'install_all'])

# The test scripts are always generated, this target allows us to
# generate the testscripts without doing a full build
env.Alias('testscripts',[])

if not IS_WINDOWS:
    generateTestScripts(env, TestGroups)



######################## Summarize our environment ###########################
def print_summary():
    print("")
    print("*** Config Summary (see config.log and lib/buildvars for details) ***")
    print("Escript/Finley revision %s"%global_revision)
    print("  Install prefix:  %s"%env['prefix'])
    print("          Python:  %s"%sysconfig.PREFIX)
    print("           boost:  %s"%env['boost_prefix'])
    print("           numpy:  YES")
    if env['usempi']:
        print("             MPI:  YES (flavour: %s)"%env['mpi'])
    else:
        print("             MPI:  DISABLED")
    if env['uselapack']:
        print("          LAPACK:  YES (flavour: %s)"%env['lapack'])
    else:
        print("          LAPACK:  DISABLED")
    d_list=[]
    e_list=[]
    for i in 'debug','openmp','boomeramg','mkl','netcdf','papi','parmetis','pyproj','silo','sympy','umfpack','visit','vsl_random':
        if env[i]: e_list.append(i)
        else: d_list.append(i)
    for i in e_list:
        print("%16s:  YES"%i)
    for i in d_list:
        print("%16s:  DISABLED"%i)
    if env['cppunit']:
        print("         CppUnit:  FOUND")
    else:
        print("         CppUnit:  NOT FOUND")
    if env['gmsh']=='m':
        print("            gmsh:  FOUND, MPI-ENABLED")
    elif env['gmsh']=='s':
        print("            gmsh:  FOUND")
    else:
        print("            gmsh:  NOT FOUND")
    if env['numpy_h']:
        print("   numpy headers:  FOUND")
    else:
        print("   numpy headers:  NOT FOUND")
    print("   vsl_random:  %s"%env['vsl_random'])
         
    if ((fatalwarning != '') and (env['werror'])):
        print("  Treating warnings as errors")
    else:
        print("  NOT treating warnings as errors")
    print("")
    for w in env['warnings']:
        print("WARNING: %s"%w)

atexit.register(print_summary)

