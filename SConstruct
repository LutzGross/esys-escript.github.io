########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

EnsureSConsVersion(0,98,1)
EnsurePythonVersion(2,5)

import sys, os, platform, re
from distutils import sysconfig
from site_init import *
import subprocess
from subprocess import PIPE, Popen

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
)

##################### Create environment and help text #######################

# Intel's compiler uses regular expressions improperly and emits a warning
# about failing to find the compilers. This warning can be safely ignored.

# PATH is needed so the compiler, linker and tools are found if they are not
# in default locations.
env = Environment(tools = ['default'], options = vars,
                  ENV = {'PATH': os.environ['PATH']})
if env['tools_names'] != 'default':
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

#################### Make sure install directories exist #####################

env['BUILD_DIR']=env['build_dir']
prefix=Dir(env['prefix']).abspath
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
    cc_flags    = "-std=c99 -fPIC -wd161 -w1 -vec-report0 -DBLOCKTIMER -DCORE_ID1"
    cc_optim    = "-O3 -ftz -IPF_ftlacc- -IPF_fma -fno-alias -ip"
    cc_debug    = "-g -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK"
    omp_flags   = "-openmp -openmp_report0"
    omp_ldflags = "-openmp -openmp_report0 -lpthread"
    fatalwarning = "-Werror"
elif cc_name[:3] == 'gcc':
    # GNU C on any system
    cc_flags     = "-pedantic -Wall -fPIC -ffast-math -Wno-unknown-pragmas -DBLOCKTIMER  -Wno-sign-compare -Wno-system-headers -Wno-long-long -Wno-strict-aliasing -finline-functions"
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

# set defaults if not otherwise specified
if env['cc_flags']    == 'default': env['cc_flags'] = cc_flags
if env['cc_optim']    == 'default': env['cc_optim'] = cc_optim
if env['cc_debug']    == 'default': env['cc_debug'] = cc_debug
if env['omp_flags']   == 'default': env['omp_flags'] = omp_flags
if env['omp_ldflags'] == 'default': env['omp_ldflags'] = omp_ldflags
if env['cc_extra']  != '': env.Append(CFLAGS = env['cc_extra'])
if env['cxx_extra'] != '': env.Append(CXXFLAGS = env['cxx_extra'])
if env['ld_extra']  != '': env.Append(LINKFLAGS = env['ld_extra'])

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
   print("OpenMP requested but no flags provided - disabling OpenMP!")
   env['openmp'] = False

if env['openmp']:
    env.Append(CCFLAGS = env['omp_flags'])
    if env['omp_ldflags'] != '': env.Append(LINKFLAGS = env['omp_ldflags'])
else:
    env['omp_flags']=''
    env['omp_ldflags']=''

# add debug/non-debug compiler flags
if env['debug']:
    env.Append(CCFLAGS = env['cc_debug'])
else:
    env.Append(CCFLAGS = env['cc_optim'])

# always add cc_flags
env.Append(CCFLAGS = env['cc_flags'])

# add system libraries
env.AppendUnique(LIBS = env['sys_libs'])


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
env.Append(CPPDEFINES=['SVN_VERSION='+global_revision])

if IS_WINDOWS:
    if not env['build_shared']:
        env.Append(CPPDEFINES = ['ESYSUTILS_STATIC_LIB'])
        env.Append(CPPDEFINES = ['PASO_STATIC_LIB'])

###################### Copy required environment vars ########################

# Windows doesn't use LD_LIBRARY_PATH but PATH instead
if IS_WINDOWS:
    LD_LIBRARY_PATH_KEY='PATH'
    env['ENV']['LD_LIBRARY_PATH']=''
else:
    LD_LIBRARY_PATH_KEY='LD_LIBRARY_PATH'

# the following env variables are exported for the unit tests

for key in 'OMP_NUM_THREADS', 'ESCRIPT_NUM_PROCS', 'ESCRIPT_NUM_NODES':
    try:
        env['ENV'][key] = os.environ[key]
    except KeyError:
        env['ENV'][key] = 1

env_export=env['env_export']
env_export.extend(['ESCRIPT_NUM_THREADS','ESCRIPT_HOSTFILE','DISPLAY','XAUTHORITY','PATH','HOME','TMPDIR','TEMP','TMP'])

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

# Create a Configure() environment to check for compilers and python
conf = Configure(env.Clone())

######## Test that the compilers work

if 'CheckCC' in dir(conf): # exists since scons 1.1.0
    if not conf.CheckCC():
        print("Cannot run C compiler '%s' (check config.log)" % (env['CC']))
        Exit(1)
    if not conf.CheckCXX():
        print("Cannot run C++ compiler '%s' (check config.log)" % (env['CXX']))
        Exit(1)
else:
    if not conf.CheckFunc('printf', language='c'):
        print("Cannot run C compiler '%s' (check config.log)" % (env['CC']))
        Exit(1)
    if not conf.CheckFunc('printf', language='c++'):
        print("Cannot run C++ compiler '%s' (check config.log)" % (env['CXX']))
        Exit(1)

if conf.CheckFunc('gethostname'):
    conf.env.Append(CPPDEFINES = ['HAVE_GETHOSTNAME'])

######## Python headers & library (required)

#First we check to see if the config file has specified
##Where to find the filae. Ideally, this should be automatic
#But we need to deal with the case where python is not in its INSTALL
#Directory
# Use the python scons is running
if env['pythoncmd']=='python':
    python_inc_path=sysconfig.get_python_inc()
    if IS_WINDOWS:
        python_lib_path=os.path.join(sysconfig.get_config_var('prefix'), 'libs')
    elif env['PLATFORM']=='darwin':
        python_lib_path=sysconfig.get_config_var('LIBPL')
    else:
        python_lib_path=sysconfig.get_config_var('LIBDIR')

    #python_libs=[sysconfig.get_config_var('LDLIBRARY')] # only on linux
    if IS_WINDOWS:
        python_libs=['python%s%s'%(sys.version_info[0], sys.version_info[1])]
    else:
        python_libs=['python'+sysconfig.get_python_version()]

#if we want to use a python other than the one scons is running
else:
    initstring='from __future__ import print_function;from distutils import sysconfig;'
    if env['pythonlibname']!='':
        python_libs=env['pythonlibname']
    else:	# work it out by calling python    
        if IS_WINDOWS:
            cmd='print("python%s%s"%(sys.version_info[0], sys.version_info[1]))'
        else:
            cmd='print("python"+sysconfig.get_python_version())'
        p=Popen([env['pythoncmd'], '-c', initstring+cmd], stdout=PIPE)
        python_libs=p.stdout.readline()
        if env['usepython3']:		# This is to convert unicode str into py2 string
            python_libs=python_libs.encode() # If scons runs on py3 then this must be rethought
        p.wait()
        python_libs=python_libs.strip()

   
    # Now we know whether we are using python3 or not
    p=Popen([env['pythoncmd'], '-c',  initstring+'print(sysconfig.get_python_inc())'], stdout=PIPE)
    python_inc_path=p.stdout.readline()
    if env['usepython3']:
         python_inc_path=python_inc_path.encode()
    p.wait()   
    python_inc_path=python_inc_path.strip()
    if IS_WINDOWS:
        cmd="os.path.join(sysconfig.get_config_var('prefix'), 'libs')"
    elif env['PLATFORM']=='darwin':
        cmd="sysconfig.get_config_var(\"LIBPL\")"
    else:
        cmd="sysconfig.get_config_var(\"LIBDIR\")"

    p=Popen([env['pythoncmd'], '-c', initstring+'print('+cmd+')'], stdout=PIPE)
    python_lib_path=p.stdout.readline()
    if env['usepython3']:
        python_lib_path=python_lib_path.decode()
    p.wait()
    python_lib_path=python_lib_path.strip()

#Check for an override from the config file.
#Ideally, this should be automatic
#But we need to deal with the case where python is not in its INSTALL
#Directory
if env['pythonlibpath']!='':
    python_lib_path=env['pythonlibpath']

if env['pythonincpath']!='':
    python_inc_path=env['pythonincpath']


if sysheaderopt == '':
    conf.env.AppendUnique(CPPPATH = [python_inc_path])
else:
    conf.env.Append(CCFLAGS = [sysheaderopt, python_inc_path])

conf.env.AppendUnique(LIBPATH = [python_lib_path])
conf.env.AppendUnique(LIBS = python_libs)
# The wrapper script needs to find the libs
conf.env.PrependENVPath(LD_LIBRARY_PATH_KEY, python_lib_path)

if not conf.CheckCHeader('Python.h'):
    print("Cannot find python include files (tried 'Python.h' in directory %s)" % (python_inc_path))
    Exit(1)
if not conf.CheckFunc('Py_Exit'):
    print("Cannot find python library method Py_Main (tried %s in directory %s)" % (python_libs, python_lib_path))
    Exit(1)

## reuse conf to check for numpy header (optional)
if env['usepython3']:
    # FIXME: This is until we can work out how to make the checks in python 3
    conf.env['numpy_h']=False
else:
    if conf.CheckCXXHeader(['Python.h','numpy/ndarrayobject.h']):
        conf.env.Append(CPPDEFINES = ['HAVE_NUMPY_H'])
        conf.env['numpy_h']=True
    else:
        conf.env['numpy_h']=False

# Commit changes to environment
env = conf.Finish()

######## boost (required)

boost_inc_path,boost_lib_path=findLibWithHeader(env, env['boost_libs'], 'boost/python.hpp', env['boost_prefix'], lang='c++')
if sysheaderopt == '':
    env.AppendUnique(CPPPATH = [boost_inc_path])
else:
    # This is required because we can't -isystem /usr/include since it breaks
    # std includes
    if os.path.normpath(boost_inc_path) == '/usr/include':
        conf.env.Append(CCFLAGS=[sysheaderopt, os.path.join(boost_inc_path,'boost')])
    else:
        env.Append(CCFLAGS=[sysheaderopt, boost_inc_path])

env.AppendUnique(LIBPATH = [boost_lib_path])
env.AppendUnique(LIBS = env['boost_libs'])
env.PrependENVPath(LD_LIBRARY_PATH_KEY, boost_lib_path)

######## numpy (required)

if env['pythoncmd']=='python':
    try:
      from numpy import identity
    except ImportError:
      print("Cannot import numpy, you need to set your PYTHONPATH and probably %s"%LD_LIBRARY_PATH_KEY)
      Exit(1)
else:
    p=subprocess.call([env['pythoncmd'],'-c','import numpy'])
    if p!=0:
      print("Cannot import numpy, you need to set your PYTHONPATH and probably %s"%LD_LIBRARY_PATH_KEY)
      Exit(1)

######## CppUnit (required for tests)

try:
    cppunit_inc_path,cppunit_lib_path=findLibWithHeader(env, env['cppunit_libs'], 'cppunit/TestFixture.h', env['cppunit_prefix'], lang='c++')
    env.AppendUnique(CPPPATH = [cppunit_inc_path])
    env.AppendUnique(LIBPATH = [cppunit_lib_path])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, cppunit_lib_path)
    env['cppunit']=True
except:
    env['cppunit']=False

######## netCDF (optional)

netcdf_inc_path=''
netcdf_lib_path=''
if env['netcdf']:
    netcdf_inc_path,netcdf_lib_path=findLibWithHeader(env, env['netcdf_libs'], 'netcdf.h', env['netcdf_prefix'], lang='c++')
    env.AppendUnique(CPPPATH = [netcdf_inc_path])
    env.AppendUnique(LIBPATH = [netcdf_lib_path])
    env.AppendUnique(LIBS = env['netcdf_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, netcdf_lib_path)
    env.Append(CPPDEFINES = ['USE_NETCDF'])

######## PAPI (optional)

papi_inc_path=''
papi_lib_path=''
if env['papi']:
    papi_inc_path,papi_lib_path=findLibWithHeader(env, env['papi_libs'], 'papi.h', env['papi_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [papi_inc_path])
    env.AppendUnique(LIBPATH = [papi_lib_path])
    env.AppendUnique(LIBS = env['papi_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, papi_lib_path)
    env.Append(CPPDEFINES = ['BLOCKPAPI'])

######## MKL (optional)

mkl_inc_path=''
mkl_lib_path=''
if env['mkl']:
    mkl_inc_path,mkl_lib_path=findLibWithHeader(env, env['mkl_libs'], 'mkl_solver.h', env['mkl_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [mkl_inc_path])
    env.AppendUnique(LIBPATH = [mkl_lib_path])
    env.AppendUnique(LIBS = env['mkl_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, mkl_lib_path)
    env.Append(CPPDEFINES = ['MKL'])

######## UMFPACK (optional)

umfpack_inc_path=''
umfpack_lib_path=''
if env['umfpack']:
    umfpack_inc_path,umfpack_lib_path=findLibWithHeader(env, env['umfpack_libs'], 'umfpack.h', env['umfpack_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [umfpack_inc_path])
    env.AppendUnique(LIBPATH = [umfpack_lib_path])
    env.AppendUnique(LIBS = env['umfpack_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, umfpack_lib_path)
    env.Append(CPPDEFINES = ['UMFPACK'])

######## LAPACK (optional)

if env['lapack']=='mkl' and not env['mkl']:
    print("mkl_lapack requires MKL!")
    Exit(1)

env['uselapack'] = env['lapack']!='none'
lapack_inc_path=''
lapack_lib_path=''
if env['uselapack']:
    header='clapack.h'
    if env['lapack']=='mkl':
        env.AppendUnique(CPPDEFINES = ['MKL_LAPACK'])
        header='mkl_lapack.h'
    lapack_inc_path,lapack_lib_path=findLibWithHeader(env, env['lapack_libs'], header, env['lapack_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [lapack_inc_path])
    env.AppendUnique(LIBPATH = [lapack_lib_path])
    env.AppendUnique(LIBS = env['lapack_libs'])
    env.Append(CPPDEFINES = ['USE_LAPACK'])

######## Silo (optional)

silo_inc_path=''
silo_lib_path=''
if env['silo']:
    silo_inc_path,silo_lib_path=findLibWithHeader(env, env['silo_libs'], 'silo.h', env['silo_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [silo_inc_path])
    env.AppendUnique(LIBPATH = [silo_lib_path])
    # Note that we do not add the libs since they are only needed for the
    # weipa library and tools.
    #env.AppendUnique(LIBS = [env['silo_libs']])

######## VSL random numbers (optional)
if env['vsl_random']:
    env.Append(CPPDEFINES = ['MKLRANDOM'])

######## VisIt (optional)

visit_inc_path=''
visit_lib_path=''
if env['visit']:
    visit_inc_path,visit_lib_path=findLibWithHeader(env, env['visit_libs'], 'VisItControlInterface_V2.h', env['visit_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [visit_inc_path])
    env.AppendUnique(LIBPATH = [visit_lib_path])

######## MPI (optional)

if env['mpi']=='no':
    env['mpi']='none'

env['usempi'] = env['mpi']!='none'
mpi_inc_path=''
mpi_lib_path=''
if env['usempi']:
    mpi_inc_path,mpi_lib_path=findLibWithHeader(env, env['mpi_libs'], 'mpi.h', env['mpi_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [mpi_inc_path])
    env.AppendUnique(LIBPATH = [mpi_lib_path])
    env.AppendUnique(LIBS = env['mpi_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, mpi_lib_path)
    env.Append(CPPDEFINES = ['ESYS_MPI', 'MPI_NO_CPPBIND', 'MPICH_IGNORE_CXX_SEEK'])
    # NetCDF 4.1 defines MPI_Comm et al. if MPI_INCLUDED is not defined!
    # On the other hand MPT and OpenMPI don't define the latter so we have to
    # do that here
    if env['netcdf'] and env['mpi'] in ['MPT','OPENMPI']:
        env.Append(CPPDEFINES = ['MPI_INCLUDED'])

######## BOOMERAMG (optional)

if env['mpi'] == 'none': env['boomeramg'] = False

boomeramg_inc_path=''
boomeramg_lib_path=''
if env['boomeramg']:
    boomeramg_inc_path,boomeramg_lib_path=findLibWithHeader(env, env['boomeramg_libs'], 'HYPRE.h', env['boomeramg_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [boomeramg_inc_path])
    env.AppendUnique(LIBPATH = [boomeramg_lib_path])
    env.AppendUnique(LIBS = env['boomeramg_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, boomeramg_lib_path)
    env.Append(CPPDEFINES = ['BOOMERAMG'])

######## ParMETIS (optional)

if not env['usempi']: env['parmetis'] = False

parmetis_inc_path=''
parmetis_lib_path=''
if env['parmetis']:
    parmetis_inc_path,parmetis_lib_path=findLibWithHeader(env, env['parmetis_libs'], 'parmetis.h', env['parmetis_prefix'], lang='c')
    env.AppendUnique(CPPPATH = [parmetis_inc_path])
    env.AppendUnique(LIBPATH = [parmetis_lib_path])
    env.AppendUnique(LIBS = env['parmetis_libs'])
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, parmetis_lib_path)
    env.Append(CPPDEFINES = ['USE_PARMETIS'])

######## gmsh (optional, for tests)

try:
    import subprocess
    p=subprocess.Popen(['gmsh', '-info'], stderr=subprocess.PIPE)
    _,e=p.communicate()
    if e.split().count("MPI"):
        env['gmsh']='m'
    else:
        env['gmsh']='s'
except OSError:
    env['gmsh']=False

######## PDFLaTeX (for documentation)
if 'PDF' in dir(env) and '.tex' in env.PDF.builder.src_suffixes(env):
    env['pdflatex']=True
else:
    env['pdflatex']=False

######################## Summarize our environment ###########################

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
for i in 'debug','openmp','netcdf','parmetis','papi','mkl','umfpack','boomeramg','silo','visit','vsl_random':
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

####################### Configure the subdirectories #########################

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

# remove obsolete file
if not env['usempi']:
    Execute(Delete(os.path.join(env['libinstall'], 'pythonMPI')))
    Execute(Delete(os.path.join(env['libinstall'], 'pythonMPIredirect')))

# Try to extract the boost version from version.hpp
boosthpp=open(os.path.join(boost_inc_path, 'boost', 'version.hpp'))
boostversion='unknown'
try:
    for line in boosthpp:
        ver=re.match(r'#define BOOST_VERSION (\d+)',line)
        if ver:
            boostversion=ver.group(1)
except StopIteration:
    pass
boosthpp.close()


buildvars=open(os.path.join(env['libinstall'], 'buildvars'), 'w')
buildvars.write("svn_revision="+str(global_revision)+"\n")
buildvars.write("prefix="+prefix+"\n")
buildvars.write("cc="+env['CC']+"\n")
buildvars.write("cxx="+env['CXX']+"\n")
if env['pythoncmd']=='python':
    buildvars.write("python="+sys.executable+"\n")
    buildvars.write("python_version="+str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2])+"\n")
else:
    buildvars.write("python="+env['pythoncmd']+"\n")
    p=Popen([env['pythoncmd'], '-c', 'from __future__ import print_function;import sys;print(str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2]))'], stdout=PIPE)
    verstring=p.stdout.readline().strip()
    p.wait()
    buildvars.write("python_version="+verstring+"\n")
buildvars.write("boost_inc_path="+boost_inc_path+"\n")
buildvars.write("boost_lib_path="+boost_lib_path+"\n")
buildvars.write("boost_version="+boostversion+"\n")
buildvars.write("debug=%d\n"%int(env['debug']))
buildvars.write("openmp=%d\n"%int(env['openmp']))
buildvars.write("mpi=%s\n"%env['mpi'])
buildvars.write("mpi_inc_path=%s\n"%mpi_inc_path)
buildvars.write("mpi_lib_path=%s\n"%mpi_lib_path)
buildvars.write("lapack=%s\n"%env['lapack'])
buildvars.write("vsl_random=%d\n"%int(env['vsl_random']))
for i in 'netcdf','parmetis','papi','mkl','umfpack','boomeramg','silo','visit':
    buildvars.write("%s=%d\n"%(i, int(env[i])))
    if env[i]:
        buildvars.write("%s_inc_path=%s\n"%(i, eval(i+'_inc_path')))
        buildvars.write("%s_lib_path=%s\n"%(i, eval(i+'_lib_path')))
buildvars.close()

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
env.Alias('run_tests', ['install_all'])
env.Alias('all_tests', ['install_all', 'run_tests', 'py_tests'])
env.Alias('build_full',['install_all','build_tests','build_py_tests'])
env.Alias('build_PasoTests','$BUILD_DIR/$PLATFORM/paso/profiling/PasoTests')

##################### Targets to build the documentation #####################

env.Alias('api_epydoc','install_all')
env.Alias('docs', ['examples_tarfile', 'examples_zipfile', 'api_epydoc', 'api_doxygen', 'user_pdf', 'install_pdf', 'cookbook_pdf'])
env.Alias('release_prep', ['docs', 'install_all'])

if not IS_WINDOWS:
    try:
        utest=open('utest.sh','w')
        utest.write(GroupTest.makeHeader(env['PLATFORM'], prefix))
        for tests in TestGroups:
            utest.write(tests.makeString())
        utest.close()
        Execute(Chmod('utest.sh', 0o755))
        print("Generated utest.sh.")
    except IOError:
        print("Error attempting to write unittests file.")
        Exit(1)

    # delete utest.sh upon cleanup
    env.Clean('target_init', 'utest.sh')

    # Make sure that the escript wrapper is in place
    if not os.path.isfile(os.path.join(env['bininstall'], 'run-escript')):
        print("Copying escript wrapper.")
        Execute(Copy(os.path.join(env['bininstall'],'run-escript'), 'bin/run-escript'))

