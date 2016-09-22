##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
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
REQUIRED_OPTS_VERSION=203

# MS Windows support, many thanks to PH
IS_WINDOWS = (os.name == 'nt')

IS_OSX = (os.uname()[0] == 'Darwin')

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
    print("We recommend that you copy the most relavent options file in the scons/template/")
    print("subdirectory and customize it to your needs.\n")
    options_file = None

############################### Build options ################################

default_prefix='/usr'
mpi_flavours=('no', 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI')
all_domains = ['dudley','finley','ripley','speckley']

#Note that scons construction vars the the following purposes:
#  CPPFLAGS -> to the preprocessor
#  CCFLAGS  -> flags for _both_ C and C++
#  CXXFLAGS -> flags for c++ _only_
#  CFLAGS   -> flags for c only

vars = Variables(options_file, ARGUMENTS)
vars.AddVariables(
  PathVariable('options_file', 'Path to options file', options_file, PathVariable.PathIsFile),
  PathVariable('prefix', 'Installation prefix', Dir('#.').abspath, PathVariable.PathIsDirCreate),
  PathVariable('build_dir', 'Top-level build directory', Dir('#/build').abspath, PathVariable.PathIsDirCreate),
  BoolVariable('verbose', 'Output full compile/link lines', False),
# Compiler/Linker options
  ('cxx', 'Path to C++ compiler', 'default'),
  ('cc_flags', 'Base (C and C++) compiler flags', 'default'),
  ('cc_optim', 'Additional (C and C++) flags for a non-debug build', 'default'),
  ('cc_debug', 'Additional (C and C++) flags for a debug build', 'default'),
  ('cxx_extra', 'Extra C++ compiler flags', ''),
  ('ld_extra', 'Extra linker flags', ''),
  ('nvcc', 'Path to CUDA compiler', 'default'),
  ('nvccflags', 'Base CUDA compiler flags', 'default'),
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
  BoolVariable('cuda', 'Enable GPU code with CUDA (requires thrust)', False),
  ('cuda_prefix', 'Prefix/Paths to NVidia CUDA installation', default_prefix),
  BoolVariable('netcdf', 'Enable netCDF file support', False),
  ('netcdf_prefix', 'Prefix/Paths of netCDF installation', default_prefix),
  ('netcdf_libs', 'netCDF libraries to link with', ['netcdf_c++', 'netcdf']),
  BoolVariable('parmetis', 'Enable ParMETIS (requires MPI)', False),
  ('parmetis_prefix', 'Prefix/Paths of ParMETIS installation', default_prefix),
  ('parmetis_libs', 'ParMETIS libraries to link with', ['parmetis', 'metis']),
  BoolVariable('mkl', 'Enable the Math Kernel Library', False),
  ('mkl_prefix', 'Prefix/Paths to MKL installation', default_prefix),
  ('mkl_libs', 'MKL libraries to link with', ['mkl_solver','mkl_em64t','guide','pthread']),
  BoolVariable('umfpack', 'Enable UMFPACK', False),
  ('umfpack_prefix', 'Prefix/Paths to UMFPACK installation', default_prefix),
  ('umfpack_libs', 'UMFPACK libraries to link with', ['umfpack']),
  BoolVariable('boomeramg', 'Enable BoomerAMG', False),
  ('boomeramg_prefix', 'Prefix/Paths to BoomerAMG installation', default_prefix),
  ('boomeramg_libs', 'BoomerAMG libraries to link with', ['boomeramg']),
  TristateVariable('lapack', 'Enable LAPACK', 'auto'),
  ('lapack_prefix', 'Prefix/Paths to LAPACK installation', default_prefix),
  ('lapack_libs', 'LAPACK libraries to link with', []),
  BoolVariable('silo', 'Enable the Silo file format in weipa', False),
  ('silo_prefix', 'Prefix/Paths to Silo installation', default_prefix),
  ('silo_libs', 'Silo libraries to link with', ['siloh5', 'hdf5']),
  BoolVariable('trilinos', 'Enable the Trilinos solvers', False),
  ('trilinos_prefix', 'Prefix/Paths to Trilinos installation', default_prefix),
  ('trilinos_libs', 'Trilinos libraries to link with', []),
  BoolVariable('visit', 'Enable the VisIt simulation interface', False),
  ('visit_prefix', 'Prefix/Paths to VisIt installation', default_prefix),
  ('visit_libs', 'VisIt libraries to link with', ['simV2']),
  ListVariable('domains', 'Which domains to build', 'all', all_domains),
  BoolVariable('paso', 'Build Paso solver library', True),
  BoolVariable('weipa', 'Build Weipa data export library', True),
# Advanced settings
  ('launcher', 'Launcher command (e.g. mpirun)', 'default'),
  ('prelaunch', 'Command to execute before launcher (e.g. mpdboot)', 'default'),
  ('postlaunch', 'Command to execute after launcher (e.g. mpdexit)', 'default'),
  #dudley_assemble_flags = -funroll-loops      to actually do something
  ('dudley_assemble_flags', 'compiler flags for some dudley optimisations', ''),
  # To enable passing function pointers through python
  BoolVariable('iknowwhatimdoing', 'Allow non-standard C', False),
  # An option for specifying the compiler tools
  ('tools_names', 'Compiler tools to use', ['default']),
  ('env_export', 'Environment variables to be passed to tools',[]),
  TristateVariable('forcelazy', 'For testing use only - set the default value for autolazy', 'auto'),
  TristateVariable('forcecollres', 'For testing use only - set the default value for force resolving collective ops', 'auto'),
  BoolVariable('build_shared', '(deprecated option, ignored)', True),
  ('sys_libs', 'Extra libraries to link with', []),
  ('escript_opts_version', 'Version of options file (do not specify on command line)'),
  ('SVN_VERSION', 'Do not use from options file', -2),
  ('pythoncmd', 'which python to compile with', sys.executable),
  ('pythonlibname', 'Name of the python library to link. (This is found automatically for python2.X.)', ''),
  ('pythonlibpath', 'Path to the python library. (You should not need to set this unless your python has moved)',''),
  ('pythonincpath','Path to python include files. (You should not need to set this unless your python has moved',''),
  BoolVariable('longindices', 'use long indices (for very large matrices)', False),
  BoolVariable('compressed_files','Enables reading from compressed binary files', True),
  ('compression_libs', 'Compression libraries to link with', ['boost_iostreams']),
  BoolVariable('papi', 'Enable PAPI', False),
  ('papi_prefix', 'Prefix/Paths to PAPI installation', default_prefix),
  ('papi_libs', 'PAPI libraries to link with', ['papi']),
  BoolVariable('papi_instrument_solver', 'Use PAPI to instrument each iteration of the solver', False),
  BoolVariable('osx_dependency_fix', 'Fix dependencies for libraries to have absolute paths (OSX)', False)
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
    env['CXX']='clang++'

if env['tools_names'] != ['default']:
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
        print("is outdated! Please update the file after reading scons/templates/README_FIRST")
        print("and setting escript_opts_version to %d.\n"%REQUIRED_OPTS_VERSION)
        Exit(1)

# Generate help text (scons -h)
Help(vars.GenerateHelpText(env))

# Check for superfluous options
if len(vars.UnknownVariables())>0:
    for k in vars.UnknownVariables():
        print("Unknown option '%s'" % k)
    Exit(1)

if env['cuda']:
    if env['nvcc'] != 'default':
        env['NVCC'] = env['nvcc']
    env.Tool('nvcc')

if 'dudley' in env['domains']:
    env['domains'].append('finley')

env['domains'] = sorted(set(env['domains']))

# create dictionary which will be populated with info for buildvars file
env['buildvars'] = {}
# create list which will be populated with warnings if there are any
env['warnings'] = []

#################### Make sure install directories exist #####################

env['BUILD_DIR'] = Dir(env['build_dir']).abspath
prefix = Dir(env['prefix']).abspath
env['buildvars']['prefix'] = prefix
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

if env['cxx'] != 'default':
    env['CXX'] = env['cxx']

# default compiler/linker options
cc_flags = '-std=c++11'
cc_optim = ''
cc_debug = ''
omp_flags = ''
omp_ldflags = ''
fatalwarning = '' # switch to turn warnings into errors
sysheaderopt = '' # how to indicate that a header is a system header

# env['CC'] might be a full path
cc_name=os.path.basename(env['CXX'])

if cc_name == 'icpc':
    # Intel compiler
    # #1478: class "std::auto_ptr<...>" was declared deprecated
    # #1875: offsetof applied to non-POD types is nonstandard (in boost)
    # removed -std=c99 because icpc doesn't like it and we aren't using c anymore
    cc_flags    = "-std=c++11 -fPIC -w2 -wd1875 -wd1478 -Wno-unknown-pragmas"
    cc_optim    = "-O3 -ftz -fno-alias -inline-level=2 -ipo -xHost"
    cc_debug    = "-g -O0 -DDOASSERT -DDOPROF -DBOUNDS_CHECK -DSLOWSHARECHECK"
    omp_flags   = "-qopenmp"
    omp_ldflags = "-qopenmp" # removing -openmp-report (which is deprecated) because the replacement outputs to a file
    fatalwarning = "-Werror"
elif cc_name[:3] == 'g++':
    # GNU C++ on any system
    # note that -ffast-math is not used because it breaks isnan(),
    # see mantis #691
    cc_flags     = "-std=c++11 -pedantic -Wall -fPIC -Wno-unknown-pragmas -Wno-sign-compare -Wno-system-headers -Wno-long-long -Wno-strict-aliasing -finline-functions"
    cc_optim     = "-O3"
    #max-vartrack-size: avoid vartrack limit being exceeded with escriptcpp.cpp
    cc_debug     = "-g3 -O0 -D_GLIBCXX_DEBUG -DDOASSERT -DDOPROF -DBOUNDS_CHECK -DSLOWSHARECHECK --param=max-vartrack-size=100000000"
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
if env['cxx_extra'] != '': env.Append(CXXFLAGS = env['cxx_extra'])
if env['ld_extra']  != '': env.Append(LINKFLAGS = env['ld_extra'])

if env['nvccflags'] != 'default':
    env['NVCCFLAGS'] = env['nvccflags']
    env['SHNVCCFLAGS'] = env['nvccflags'] + ' -shared'

if env['longindices']:
    env.Append(CPPDEFINES = ['ESYS_INDEXTYPE_LONG'])

# set up the autolazy values
if env['forcelazy'] == 1:
    env.Append(CPPDEFINES=['FAUTOLAZYON'])
elif env['forcelazy'] == 0:
    env.Append(CPPDEFINES=['FAUTOLAZYOFF'])

# set up the collective resolve values
if env['forcecollres'] == 1:
    env.Append(CPPDEFINES=['FRESCOLLECTON'])
elif env['forcecollres'] == 0:
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

env['IS_WINDOWS']=IS_WINDOWS
env['IS_OSX']=IS_OSX

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
        env['ENV'][key] = '1'

env_export=env['env_export']
env_export.extend(['ESCRIPT_NUM_THREADS','ESCRIPT_HOSTFILE','DISPLAY','XAUTHORITY','PATH','HOME','KMP_MONITOR_STACKSIZE','TMPDIR','TEMP','TMP','LD_PRELOAD'])

for key in set(env_export):
    try:
        env['ENV'][key] = os.environ[key]
    except KeyError:
        pass

for key in os.environ.keys():
    if key.startswith("SLURM_"):
        env['ENV'][key] = os.environ[key]

try:
    env.PrependENVPath(LD_LIBRARY_PATH_KEY, os.environ[LD_LIBRARY_PATH_KEY])
except KeyError:
    pass

if IS_OSX:
  try:
    env.PrependENVPath('DYLD_LIBRARY_PATH', os.environ['DYLD_LIBRARY_PATH'])
  except KeyError:
    pass

try:
    env['ENV']['PYTHONPATH'] = os.environ['PYTHONPATH']
except KeyError:
    pass

######################## Add some custom builders ############################

# Takes care of prefix and suffix for Python modules:
def build_python_module(env, target, source):
    return env.SharedLibrary(target, source, SHLIBPREFIX='', SHLIBSUFFIX='.so')
env.AddMethod(build_python_module, "PythonModule")

if env['pythoncmd']=='python':
    py_builder = Builder(action = build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
else:
    py_builder = Builder(action = env['pythoncmd']+" scripts/py_comp.py $SOURCE $TARGET", suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

runUnitTest_builder = Builder(action = runUnitTest, suffix = '.passed', src_suffix=env['PROGSUFFIX'], single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

runPyExample_builder = Builder(action = runPyExample, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyExample' : runPyExample_builder});

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

######## NVCC version (optional)
if env['cuda'] and 'ripley' in env['domains']:
    env=checkCudaVersion(env)
    env=checkCUDA(env)

######## optional python modules (sympy, pyproj)
env=checkOptionalModules(env)

######## optional dependencies (netCDF, PAPI, MKL, UMFPACK, Lapack, Silo, ...)
env=checkOptionalLibraries(env)

######## PDFLaTeX (for documentation)
env=checkPDFLatex(env)

# set defaults for launchers if not otherwise specified
if env['prelaunch'] == 'default':
    if env['mpi'] == 'INTELMPI' and env['openmp']:
        env['prelaunch'] = "export I_MPI_PIN_DOMAIN=omp"
    elif env['mpi'] == 'OPENMPI':
        # transform comma-separated list to '-x a -x b -x c ...'
        env['prelaunch'] = "EE=$(echo -x %e|sed -e 's/,/ -x /g')"
    elif env['mpi'] == 'MPT':
        env['prelaunch'] = "export MPI_NUM_MEMORY_REGIONS=0"
    elif env['mpi'] == 'MPICH2':
        env['prelaunch'] = "mpdboot -n %n -r ssh -f %f"
    else:
        env['prelaunch'] = ""

if env['launcher'] == 'default':
    if env['mpi'] == 'INTELMPI':
        env['launcher'] = "mpirun -hostfile %f -n %N -ppn %p %b"
    elif env['mpi'] == 'OPENMPI':
        # default to OpenMPI version 1.10 or higher
        env['launcher'] = "mpirun ${AGENTOVERRIDE} --gmca mpi_warn_on_fork 0 ${EE} --host %h --map-by node:pe=%t -bind-to core -np %N %b"
        if 'orte_version' in env:
            major,minor,point = [int(i) for i in env['orte_version'].split('.')]
            if major == 1 and minor < 10:
                env['launcher'] = "mpirun ${AGENTOVERRIDE} --gmca mpi_warn_on_fork 0 ${EE} --host %h --cpus-per-rank %t -np %N %b"
    elif env['mpi'] == 'MPT':
        env['launcher'] = "mpirun %h -np %p %b"
    elif env['mpi'] == 'MPICH':
        env['launcher'] = "mpirun -machinefile %f -np %N %b"
    elif env['mpi'] == 'MPICH2':
        env['launcher'] = "mpiexec -genvlist %e -np %N %b"
    else:
        env['launcher'] = "%b"

if env['postlaunch'] == 'default':
    if env['mpi'] == 'MPICH2':
        env['postlaunch'] = "mpdallexit"
    else:
        env['postlaunch'] = ""

# dependency sanity checks

if len(env['domains']) == 0:
   env['warnings'].append("No domains have been built, escript will not be very useful!")

# keep some of our install paths first in the list for the unit tests
env.PrependENVPath(LD_LIBRARY_PATH_KEY, env['libinstall'])
env.PrependENVPath('PYTHONPATH', prefix)
env['ENV']['ESCRIPT_ROOT'] = prefix

if not env['verbose']:
    env['CXXCOMSTR'] = "Compiling $TARGET"
    env['SHCXXCOMSTR'] = "Compiling $TARGET"
    env['ARCOMSTR'] = "Linking $TARGET"
    env['LINKCOMSTR'] = "Linking $TARGET"
    env['SHLINKCOMSTR'] = "Linking $TARGET"
    env['PDFLATEXCOMSTR'] = "Building $TARGET from LaTeX input $SOURCES"
    env['BIBTEXCOMSTR'] = "Generating bibliography $TARGET"
    env['MAKEINDEXCOMSTR'] = "Generating index $TARGET"
    env['PDFLATEXCOMSTR'] = "Building $TARGET from LaTeX input $SOURCES"
    #Progress(['Checking -\r', 'Checking \\\r', 'Checking |\r', 'Checking /\r'], interval=17)

########################### Configure the targets ############################

from grouptest import GroupTest
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

target_init = env.Command(os.path.join(env['pyinstall'],'__init__.py'), None, Touch('$TARGET'))
env.Alias('target_init', [target_init])

# escript can't be turned off
build_all_list = ['build_escript']
install_all_list = ['target_init', 'install_escript']

if env['usempi']:
    build_all_list += ['build_pythonMPI', 'build_overlord']
    install_all_list += ['install_pythonMPI', 'install_overlord']

env['buildvars']['paso'] = int(env['paso'])
if env['paso']:
    env.Append(CPPDEFINES = ['ESYS_HAVE_PASO'])
    build_all_list += ['build_paso']
    install_all_list += ['install_paso']

env['buildvars']['paso'] = int(env['trilinos'])
if env['trilinos']:
    build_all_list += ['build_trilinoswrap']
    install_all_list += ['install_trilinoswrap']

env['buildvars']['domains'] = ','.join(env['domains'])
for domain in env['domains']:
    env.Append(CPPDEFINES = ['ESYS_HAVE_'+domain.upper()])
    build_all_list += ['build_%s'%domain]
    install_all_list += ['install_%s'%domain]

env['buildvars']['weipa'] = int(env['weipa'])
if env['weipa']:
    env.Append(CPPDEFINES = ['ESYS_HAVE_WEIPA'])
    build_all_list += ['build_weipa']
    install_all_list += ['install_weipa']
    if 'finley' in env['domains'] or 'dudley' in env['domains']:
        build_all_list += ['build_escriptreader']
        install_all_list += ['install_escriptreader']

variant='$BUILD_DIR/$PLATFORM/'
env.SConscript('escriptcore/SConscript', variant_dir=variant+'escriptcore', duplicate=0)
env.SConscript('escript/py_src/SConscript', variant_dir=variant+'escript', duplicate=0)
env.SConscript('pythonMPI/src/SConscript', variant_dir=variant+'pythonMPI', duplicate=0)
env.SConscript('tools/overlord/SConscript', variant_dir=variant+'tools/overlord', duplicate=0)
env.SConscript('paso/SConscript', variant_dir=variant+'paso', duplicate=0)
env.SConscript('trilinoswrap/SConscript', variant_dir=variant+'trilinoswrap', duplicate=0)
env.SConscript('cusplibrary/SConscript')
env.SConscript('dudley/SConscript', variant_dir=variant+'dudley', duplicate=0)
env.SConscript('finley/SConscript', variant_dir=variant+'finley', duplicate=0)
env.SConscript('ripley/SConscript', variant_dir=variant+'ripley', duplicate=0)
env.SConscript('speckley/SConscript', variant_dir=variant+'speckley', duplicate=0)
env.SConscript('weipa/SConscript', variant_dir=variant+'weipa', duplicate=0)
env.SConscript(dirs = ['downunder/py_src'], variant_dir=variant+'downunder', duplicate=0)
env.SConscript(dirs = ['modellib/py_src'], variant_dir=variant+'modellib', duplicate=0)
env.SConscript(dirs = ['pycad/py_src'], variant_dir=variant+'pycad', duplicate=0)
env.SConscript('tools/escriptconvert/SConscript', variant_dir=variant+'tools/escriptconvert', duplicate=0)
env.SConscript('doc/SConscript', variant_dir=variant+'doc', duplicate=0)

env.Alias('build', build_all_list)

install_all_list += ['install_downunder_py']
install_all_list += ['install_modellib_py']
install_all_list += ['install_pycad_py']
install_all_list += [env.Install(Dir('scripts',env['build_dir']), os.path.join('scripts', 'release_sanity.py'))]

if env['osx_dependency_fix']:
    print("Require dependency fix")
    install_all=env.Command('install', install_all_list, 'scripts/moveall.sh')
else:
    install_all=env.Alias('install', install_all_list)

sanity=env.Alias('sanity', env.Command('dummy','',os.path.join(env['prefix'], 'bin', 'run-escript')+' '+os.path.join(env['build_dir'],'scripts', 'release_sanity.py')))
env.Depends('dummy', install_all)
if env['usempi']:
   env.Depends('dummy', ['install_pythonMPI'])

# if all domains are built:
if env['domains'] == all_domains:
       env.AlwaysBuild('sanity')
       env.Default('sanity')
else:
    env.Default('install')

################## Targets to build and run the test suite ###################

if not env['cppunit']:
    test_msg = env.Command('.dummy.', None, '@echo "Cannot run C++ unit tests, CppUnit not found!";exit 1')
    env.Alias('run_tests', test_msg)
    env.Alias('build_tests', '')
env.Alias('run_tests', ['install'])
env.Alias('all_tests', ['install', 'run_tests', 'py_tests'])
env.Alias('build_full',['install','build_tests','build_py_tests'])
Requires('py_tests', 'install')

##################### Targets to build the documentation #####################

env.Alias('pdfdocs',['user_pdf', 'install_pdf', 'cookbook_pdf', 'inversion_pdf'])
env.Alias('basedocs', ['pdfdocs','examples_tarfile', 'examples_zipfile', 'api_doxygen'])
env.Alias('docs', ['basedocs', 'sphinxdoc'])
env.Alias('release_prep', ['docs', 'install'])
env.Alias('release_prep_old', ['basedocs', 'api_epydoc', 'install'])

# The test scripts are always generated, this target allows us to
# generate the testscripts without doing a full build
env.Alias('testscripts',[])

if not IS_WINDOWS:
    generateTestScripts(env, TestGroups)

######################## Populate the buildvars file #########################

write_buildvars(env)
# delete buildvars upon cleanup - target_init is default so use it
env.Clean('target_init', File('buildvars', env['libinstall']))

write_launcher(env)

# remove obsolete files
if not env['usempi']:
    Execute(Delete(File(['pythonMPI','pythonMPIredirect'], env['libinstall'])))
    Execute(Delete(File('escript-overlord', env['bininstall'])))

######################## Summarize our environment ###########################
def print_summary():
    d_list=[]
    print("")
    print("*** Config Summary (see config.log and <prefix>/lib/buildvars for details) ***")
    print("Escript revision %s"%global_revision)
    print("  Install prefix:  %s"%env['prefix'])
    print("          Python:  %s (Version %s)"%(env['pythoncmd'],env['python_version']))
    print("           boost:  %s (Version %s)"%(env['boost_prefix'],env['boost_version']))
    if env['numpy_h']:
        print("           numpy:  YES (with headers)")
    else:
        print("           numpy:  YES (without headers)")
    if env['usempi']:
        if 'orte_version' in env:
            print("             MPI:  %s (Version %s)"%(env['mpi'], env['orte_version']))
        else:
            print("             MPI:  YES (flavour: %s)"%env['mpi'])
    else:
        d_list.append('mpi')
    if env['parmetis']:
        print("        ParMETIS:  %s (Version %s)"%(env['parmetis_prefix'],env['parmetis_version']))
    else:
        d_list.append('parmetis')
    if env['uselapack']:
        print("          LAPACK:  YES (flavour: %s)"%env['lapack'])
    else:
        d_list.append('lapack')
    if env['cuda']:
        print("            CUDA:  YES (nvcc: %s)"%env['nvcc_version'])
    else:
        d_list.append('cuda')
    if env['gmshpy']:
        gmshpy=" + python module"
    else:
        gmshpy=""
    if env['gmsh']=='m':
        print("            gmsh:  YES, MPI-ENABLED"+gmshpy)
    elif env['gmsh']=='s':
        print("            gmsh:  YES"+gmshpy)
    else:
        if env['gmshpy']:
            print("            gmsh:  python module only")
        else:
            d_list.append('gmsh')
    if env['compressed_files']:
        print("            gzip:  YES")
    else:
        d_list.append('gzip')

    solvers = []
    direct = []
    if env['paso']:
        solvers.append('paso')
        if env['mkl']:
            direct.append('mkl')
        if env['umfpack']:
            direct.append('umfpack')
    else:
        d_list.append('paso')
    if env['trilinos']:
        solvers.append('trilinos')
        direct.append('trilinos')
    else:
        d_list.append('trilinos')

    print("  Solver library:  %s"%(", ".join(solvers)))
    if len(direct) > 0:
        print("   Direct solver:  YES (%s)"%(", ".join(direct)))
    else:
        print("   Direct solver:  NONE")
    print("         domains:  %s"%(", ".join(env['domains'])))

    e_list=[]
    for i in 'weipa','debug','openmp','boomeramg','cppunit','gdal','mkl',\
             'netcdf','papi','pyproj','scipy','silo','sympy','umfpack','visit':
        if env[i]: e_list.append(i)
        else: d_list.append(i)

    d_list += set(all_domains).difference(env['domains'])
    for i in e_list:
        print("%16s:  YES"%i)
    print("\n  DISABLED features: %s"%(" ".join(sorted(d_list))))

    if ((fatalwarning != '') and (env['werror'])):
        print("  Treating warnings as errors")
    else:
        print("  NOT treating warnings as errors")
    print("")
    for w in env['warnings']:
        print("WARNING: %s"%w)
    if len(GetBuildFailures()):
        print("\nERROR: build stopped due to errors\n")
    else:
        print("\nSUCCESS: build complete\n")

atexit.register(print_summary)

