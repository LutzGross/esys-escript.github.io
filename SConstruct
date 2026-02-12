##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

EnsureSConsVersion(0,98,1)
EnsurePythonVersion(3,1)

import atexit, sys, os, platform, re, shutil
#from distutils import sysconfig
from dependencies import *
from site_init import *


print(sys.version)

# Version number to check for in options file. Increment when new features are
# added or existing options changed.
REQUIRED_OPTS_VERSION=203

# MS Windows support, many thanks to PH
IS_WINDOWS = (os.name == 'nt')

if IS_WINDOWS:
    IS_OSX = False
else:
    IS_OSX = (os.uname()[0] == 'Darwin')

########################## Determine options file ############################
# 1. command line
# 2. scons/<hostname>_options.py
# 3. name as part of a cluster
options_file=ARGUMENTS.get('options_file', None)
if not options_file:
    ext_dir = os.path.join(os.getcwd(), 'scons')
    hostname = platform.node().split('.')[0]
    mangledhostname = re.sub('[^0-9a-zA-Z]', '_', hostname)
    options_file = os.path.join(ext_dir, mangledhostname+'_options.py')

if not os.path.isfile(options_file):
    print("\nWARNING:\nOptions file %s" % options_file)
    print("not found! Default options will be used which is most likely suboptimal.")
    print("We recommend that you copy the most relavent options file in the scons/template/")
    print("subdirectory and customize it to your needs.\n")
    options_file = None

############################### Build options ################################

default_prefix='/usr'
mpi_flavours=('auto', 'none', 'MPT', 'MPICH', 'MPICH2', 'OPENMPI', 'INTELMPI')
# Dynamically determine available domains based on which directories exist
# This allows release tarballs to exclude optional domains like oxley
_possible_domains = ('finley','oxley','ripley','speckley')
all_domains = tuple(d for d in _possible_domains if os.path.isdir(d))
version_info=['0.0','5','6']
build_trilinos_flavours = ( "check",      # check if rebuild needed (default - fast)
                            "make",       # same as check (for backward compatibility)
                            "always",     # always rebuild (slow - ignores timestamps)
                            "never"       # never build (use existing installation)
                            )

#Note that scons construction vars the following purposes:
#  CPPFLAGS -> to the preprocessor
#  CCFLAGS  -> flags for _both_ C and C++
#  CXXFLAGS -> flags for c++ _only_
#  CFLAGS   -> flags for c only

vars = Variables(options_file, ARGUMENTS)
vars.AddVariables(
  PathVariable('options_file', 'Path to options file', options_file, PathVariable.PathIsFile),
  PathVariable('prefix', 'Installation prefix (may not contain space)', Dir('#.').abspath, PathVariable.PathIsDirCreate),
  PathVariable('PREFIX', 'Installation prefix  (may not contain space)', Dir('#.').abspath, PathVariable.PathIsDirCreate),
  PathVariable('build_dir', 'Top-level build directory', Dir('#/build').abspath, PathVariable.PathIsDirCreate),
  BoolVariable('verbose', 'Output full compile/link lines', False),
  BoolVariable('skip_link_checks', 'Skip library link checks during configuration (useful for Homebrew/macOS)', False),
# Compiler/Linker options
  ('cxx', 'Path to C++ compiler', 'default'),
  ('cc', 'Path to C compiler', 'default'),
  ('cc_flags', 'Base (C and C++) compiler flags', 'default'),
    ('cxx_flags', 'Base C++ compiler flags (if [] set to cc_flags)', []),
    ('cc_optim', 'Additional (C and C++) flags for a non-debug build', 'default'),
  ('cc_debug', 'Additional (C and C++) flags for a debug build', 'default'),
  ('cxx_extra', 'Extra C++ compiler flags', [] ),
    ('cc_extra', 'Extra C compiler flags', [] ),
    ('ld_flags', 'linker flags', 'default'),
    ('ld_extra', 'Extra linker flags', [] ),
  BoolVariable('werror', 'Treat compiler warnings as errors', True),
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

  BoolVariable('sympy', 'Enable sympy symbolic module support', False),
  BoolVariable('hdf5', 'Enable hdf5, if available', True),
  ('hdf5_prefix', 'Prefix/Paths of hdf5 installation', default_prefix),
  ('hdf5_libs', 'HDF5 libraries to link with', 'DEFAULT'),
    BoolVariable('zlib', 'Enable zLib', False),
    ('zlib_prefix', 'Prefix/Paths to zlib installation', default_prefix),
    ('zlib_libs', 'zlib libraries to link with', ['zlib']),
    BoolVariable('metis', 'Enable METIS', False),
  ('metis_prefix', 'Prefix/Paths of METIS installation', default_prefix),
  ('metis_libs', 'METIS libraries to link with', ['metis']),
    BoolVariable('parmetis', 'Enable ParMETIS (requires MPI)', False),
  ('parmetis_prefix', 'Prefix/Paths of ParMETIS installation', default_prefix),
  ('parmetis_libs', 'ParMETIS libraries to link with', ['parmetis']),
    BoolVariable('scotch', 'Enable Scotch/PT-Scotch graph partitioning', False),
  ('scotch_prefix', 'Prefix/Paths of Scotch installation', default_prefix),
  ('scotch_libs', 'Scotch libraries to link with', ['ptscotch', 'ptscotcherr', 'scotch', 'scotcherr']),
  BoolVariable('mkl', 'Enable the Math Kernel Library', False),
  ('mkl_prefix', 'Prefix/Paths to MKL installation', default_prefix),
  ('mkl_libs', 'MKL libraries to link with', ['mkl_solver','mkl_em64t','guide','pthread']),
  BoolVariable('umfpack', 'Enable UMFPACK', False),
  ('umfpack_prefix', 'Prefix/Paths to UMFPACK installation', default_prefix),
  ('umfpack_libs', 'UMFPACK libraries to link with', ['umfpack']),
  BoolVariable('mumps_seq', 'Enable sequential MUMPS (works with MPI builds)', False),
  ('mumps_seq_prefix', 'Prefix/Paths to sequential MUMPS installation', default_prefix),
  ('mumps_seq_libs', 'Sequential MUMPS libraries to link with', ['mumps_common','pord','dmumps','zmumps',
    'mpiseq','lapack','metis','scotch','esmumps','gfortran']),
  TristateVariable('lapack', 'Enable LAPACK', 'auto'),
  ('lapack_prefix', 'Prefix/Paths to LAPACK installation', default_prefix),
  ('lapack_libs', 'LAPACK libraries to link with', []),
  BoolVariable('silo', 'Enable the Silo file format in weipa', False),
  ('silo_prefix', 'Prefix/Paths to Silo installation', default_prefix),
  ('silo_libs', 'Silo libraries to link with', ['siloh5', 'hdf5']),
  ('netcdf', 'Enable NetCDF (4 for NetCDF4, False to disable)', False),
  ('netcdf_prefix', 'Prefix/Paths to NetCDF installation', default_prefix),
  ('netcdf_libs', 'NetCDF libraries to link with', ['netcdf_c++4', 'netcdf']),
  BoolVariable('trilinos', 'Enable the Trilinos solvers (overwritten when trilinos is built)', False),
  EnumVariable('build_trilinos', 'Instructs scons to build the trilinos library.', "make", allowed_values = build_trilinos_flavours),
  ('trilinos_prefix', 'Prefix/Paths to Trilinos installation (need to be set  if build_trilinos = False).', default_prefix),
  ('trilinos_libs', 'Trilinos libraries to link with', []),
  PathVariable('trilinos_src', 'Top-level source directory for trilinos.', Dir('trilinos_source_17-0-0').abspath, PathVariable.PathIsDir),
  PathVariable('trilinos_build', 'Top-level build directory for trilinos.', Dir('#/build_trilinos').abspath, PathVariable.PathIsDirCreate),
  PathVariable('trilinos_install', 'Top-level install directory for trilinos when built', Dir('#/esys.trilinos').abspath, PathVariable.PathIsDirCreate),
  #('trilinos_install', 'path to install trilinos libs, default is <prefix>/lib/esys', 'default'),
  ('trilinos_make_sh', 'path to a shell script to run trilinos make.', 'default'),
  BoolVariable('visit', 'Enable the VisIt simulation interface', False),
  ('visit_prefix', 'Prefix/Paths to VisIt installation', default_prefix),
  ('visit_libs', 'VisIt libraries to link with', ['simV2']),
  #ListVariable('domains', 'Which domains to build', 'all', all_domains),
  ('domains', 'which domains to build', all_domains),
  BoolVariable('paso', 'Build Paso solver library', True),
  BoolVariable('weipa', 'Build Weipa data export library', True),
  ('mathjax_path', 'Path to MathJax.js file', 'default'),
# Advanced settings
  ('launcher', 'Launcher command (e.g. mpirun)', 'default'),
  ('prelaunch', 'Command to execute before launcher (e.g. mpdboot)', 'default'),
  ('postlaunch', 'Command to execute after launcher (e.g. mpdexit)', 'default'),
  # To enable passing function pointers through python
  BoolVariable('iknowwhatimdoing', 'Allow non-standard C', False),
  # An option for specifying the compiler tools
  ('tools_names', 'Compiler tools to use', ['default']),
  ('env_export', 'Environment variables to be passed to tools',[]),
  TristateVariable('forcelazy', 'For testing use only - set the default value for autolazy', 'auto'),
  TristateVariable('forcecollres', 'For testing use only - set the default value for force resolving collective ops', 'auto'),
  BoolVariable('build_shared', '(deprecated option, ignored)', True),
  ('sys_libs', 'Extra libraries to link with', []),
  ('sysheaderopt', 'system header options', ''),
  ('escript_opts_version', 'Version of options file (do not specify on command line)'),
  ('SVN_VERSION', 'Do not use from options file', -2),
  ('pythoncmd', 'which python to compile with', sys.executable),
  ('pythonlibname', 'Name of the python library to link. (This is found automatically for python3.X.)', ''),
  ('pythonlibpath', 'Path to the python library. (You should not need to set this unless your python has moved)',''),
  ('pythonincpath', 'Path to python include files. (You should not need to set this unless your python has moved',''),
  BoolVariable('longindices', 'use long indices (for very large matrices)', False),
  BoolVariable('compressed_files', 'Enables reading from compressed binary files', True),
  ('compression_libs', 'Compression libraries to link with', ['boost_iostreams']),
  BoolVariable('disable_boost_numpy', 'Do not build using boost_numpy, even if it is available', False),
  BoolVariable('osx_dependency_fix', 'Fix dependencies for libraries to have absolute paths (OSX)', False),
  BoolVariable('stdlocationisprefix', 'Set the prefix as escript root in the launcher', False),
  BoolVariable('mpi_no_host', 'Do not specify --host in run-escript launcher (only OPENMPI)', False),
  BoolVariable('insane', 'Instructs scons to not run a sanity check after compilation.', False),
  EnumVariable('version_information', 'Instructs scons to create symlinks to the library files','0.0',allowed_values=version_info),
  BoolVariable('mpi4py', 'Compile with mpi4py.', False),
  BoolVariable('p4est', 'Compile with p4est.', True),
  ('trilinos_LO', 'Manually specify the LO used by Trilinos.', ''),
  ('trilinos_GO', 'Manually specify the GO used by Trilinos.', '')
)


##################### Create environment and help text #######################

# Intel's compiler uses regular expressions improperly and emits a warning
# about failing to find the compilers. This warning can be safely ignored.

# PATH is needed so the compiler, linker and tools are found if they are not
# in default locations.
#env = Environment(tools = ['default'] + [ 'applelink' ], options = vars,


env = Environment(tools = ['default'] , options = vars,
                  ENV = {'PATH': os.environ['PATH']})


# check options file version:
if options_file:
    opts_valid=False
    if 'escript_opts_version' in env.Dictionary() and \
        int(env['escript_opts_version']) >= REQUIRED_OPTS_VERSION:
            opts_valid=True
    if opts_valid:
        print("Using options in %s." % options_file)
    else:
        print("\nOptions file %s" % options_file)
        print("is outdated! Please update the file after reading scons/templates/options.md")
        print("and setting escript_opts_version to %d.\n"%REQUIRED_OPTS_VERSION)
        Exit(1)

############### Migrate deprecated option values ####################
# Migrate mpi='no' to 'none'
if env['mpi'] == 'no':
    print("WARNING: mpi='no' is deprecated, using 'none' instead. Please update your options file.")
    env['mpi'] = 'none'

# Migrate build_trilinos string-boolean values
if env['build_trilinos'] == "True":
    print("WARNING: build_trilinos='True' is deprecated, using 'make' instead. Please update your options file.")
    env['build_trilinos'] = 'make'
elif env['build_trilinos'] == "False":
    print("WARNING: build_trilinos='False' is deprecated, using 'never' instead. Please update your options file.")
    env['build_trilinos'] = 'never'
elif env['build_trilinos'] is True:
    print("WARNING: build_trilinos=True (boolean) is deprecated, using 'make' instead. Please update your options file.")
    env['build_trilinos'] = 'make'
elif env['build_trilinos'] is False:
    print("WARNING: build_trilinos=False (boolean) is deprecated, using 'never' instead. Please update your options file.")
    env['build_trilinos'] = 'never'

# Migrate integer boolean values (0/1) to True/False for consistency
# SCons accepts these but True/False is clearer
bool_options = ['openmp', 'paso', 'weipa', 'sympy', 'hdf5', 'debug', 'verbose',
                'werror', 'trilinos', 'umfpack', 'mkl', 'mumps_seq', 'silo', 'visit',
                'metis', 'parmetis', 'scotch', 'zlib', 'p4est', 'mpi4py', 'longindices',
                'compressed_files', 'disable_boost_numpy', 'osx_dependency_fix',
                'stdlocationisprefix', 'mpi_no_host', 'insane', 'iknowwhatimdoing']
for opt in bool_options:
    if opt in env and env[opt] == 1:
        env[opt] = True
    elif opt in env and env[opt] == 0:
        env[opt] = False

env['IS_WINDOWS']=IS_WINDOWS
env['IS_OSX']=IS_OSX
################# Fill in compiler options if not set above ##################
if env['cxx'] != 'default':
    env['CXX'] = env['cxx']
if env['cc'] != 'default':
    env['CC'] = env['cc']
if ( env['cc'] == 'default' and env['cxx'] != 'default') :
    env['CC'] = env['cxx']

# Auto-detect MPI flavour from mpi4py if mpi='auto'
if env['mpi'] == 'auto':
    if env['mpi4py']:
        try:
            from mpi4py import MPI
            mpi_version = MPI.Get_library_version()
            if 'Open MPI' in mpi_version:
                env['mpi'] = 'OPENMPI'
                print("Auto-detected MPI flavour from mpi4py: OPENMPI")
            elif 'MPICH' in mpi_version or 'MVAPICH' in mpi_version:
                env['mpi'] = 'MPICH'
                print("Auto-detected MPI flavour from mpi4py: MPICH")
            elif 'Intel(R) MPI' in mpi_version:
                env['mpi'] = 'INTELMPI'
                print("Auto-detected MPI flavour from mpi4py: INTELMPI")
            else:
                print("ERROR: Could not auto-detect MPI flavour from mpi4py.")
                print(f"       MPI library version: {mpi_version}")
                print("       Please set mpi='OPENMPI', 'MPICH', or 'INTELMPI' explicitly.")
                Exit(1)
        except ImportError:
            print("ERROR: mpi='auto' with mpi4py=True but mpi4py is not installed.")
            print("       Please install mpi4py or set mpi flavour explicitly.")
            Exit(1)
    else:
        # mpi4py not enabled, disable MPI
        env['mpi'] = 'none'
        print("MPI flavour set to 'none' (mpi4py not enabled)")

# Verify MPI flavour matches mpi4py if both are enabled
if env['mpi4py'] and env['mpi'] != 'none':
    try:
        from mpi4py import MPI
        mpi_version = MPI.Get_library_version()
        # Check for compatibility
        if env['mpi'] == 'OPENMPI' and 'Open MPI' not in mpi_version:
            print(f"ERROR: mpi='OPENMPI' but mpi4py is built against: {mpi_version}")
            print("       MPI flavour must match mpi4py's MPI implementation.")
            Exit(1)
        elif env['mpi'] in ('MPICH', 'MPICH2') and 'MPICH' not in mpi_version and 'MVAPICH' not in mpi_version:
            print(f"ERROR: mpi='{env['mpi']}' but mpi4py is built against: {mpi_version}")
            print("       MPI flavour must match mpi4py's MPI implementation.")
            Exit(1)
        elif env['mpi'] == 'INTELMPI' and 'Intel(R) MPI' not in mpi_version:
            print(f"ERROR: mpi='INTELMPI' but mpi4py is built against: {mpi_version}")
            print("       MPI flavour must match mpi4py's MPI implementation.")
            Exit(1)
    except ImportError:
        print("ERROR: mpi4py=True but mpi4py is not installed.")
        print("       Please install mpi4py or set mpi4py=False.")
        Exit(1)

if env['mpi'] == 'OPENMPI':
    env['CXX'] = 'mpic++'
    env['CC'] = 'mpicc'

if env['mpi4py'] and env['mpi'] == 'none':
    print("ERROR: mpi4py is enabled but mpi='none'.")
    print("       Set mpi='auto' to auto-detect or specify mpi flavour explicitly.")
    Exit(1)
# set the vars for clang
def mkclang(env):
        env['CXX']='clang'

if env['tools_names'] != ['default']:
    zz=env['tools_names']
    if 'clang' in zz:
        zz.remove('clang')
        zz.insert(0, mkclang)
        env = Environment(tools = ['default'] + env['tools_names'], options = vars,
                      ENV = {'PATH' : os.environ['PATH']})
    else:
        env = Environment(tools = ['default'] + env['tools_names'], options = vars,
                      ENV = {'PATH' : os.environ['PATH']} )


# Generate help text (scons -h)
Help(vars.GenerateHelpText(env))
# Check for superfluous options
if len(vars.UnknownVariables())>0:
    APPROVED_VARIABLES = [ 'os', 'subprocess']
    DEL= ""
    for k in vars.UnknownVariables():
        if not k in APPROVED_VARIABLES:
            print("Unknown option variable '%s'" % k)
            DEL+=k+", "
    if len(DEL)>0:
        print("Most likely these variables are set in the options file. Add the following statement to delete them:")
        print("\n  del "+DEL[:-2]+"\n")
        Exit(1)

env['domains'] = sorted(set(env['domains']))
#===== First thing we do is to install Trilinos if requested to do so:
#
################ If requested, build & install Trilinos ####################
# Manually change the trilinos ordinals (if necessary)
if env['trilinos_LO'] != '':
    env.Append(CPPDEFINES=['MANUALLY_SET_LO'])
    print("Manually setting the Trilinos Local Ordinate...")
    if env['trilinos_LO'] == 'int':
        env.Append(CPPDEFINES=['SET_LO_INT'])
    elif env['trilinos_LO'] == 'long':
        env.Append(CPPDEFINES=['SET_LO_LONG'])
    elif env['trilinos_LO'] == 'long long':
        env.Append(CPPDEFINES=['SET_LO_LONG_LONG'])
    elif env['trilinos_LO'] == 'complex double':
        env.Append(CPPDEFINES=['SET_LO_COMPLEX_DOUBLE'])
    elif env['trilinos_LO'] == 'real_t':
        env.Append(CPPDEFINES=['SET_LO_REALT'])
    elif env['trilinos_LO'] == 'cplx_t':
        env.Append(CPPDEFINES=['SET_LO_CPLXT'])
if env['trilinos_GO'] != '':
    env.Append(CPPDEFINES=['MANUALLY_SET_GO'])
    print("Manually setting the Trilinos Global Ordinate...")
    if env['trilinos_GO'] == 'int':
        env.Append(CPPDEFINES=['SET_GO_INT'])
    elif env['trilinos_GO'] == 'long':
        env.Append(CPPDEFINES=['SET_GO_LONG'])
    elif env['trilinos_GO'] == 'long long':
        env.Append(CPPDEFINES=['SET_GO_LONG_LONG'])
    elif env['trilinos_GO'] == 'complex double':
        env.Append(CPPDEFINES=['SET_GO_COMPLEX_DOUBLE'])
    elif env['trilinos_GO'] == 'real_t':
        env.Append(CPPDEFINES=['SET_GO_REALT'])
    elif env['trilinos_GO'] == 'cplx_t':
        env.Append(CPPDEFINES=['SET_GO_CPLXT'])

# create dictionary which will be populated with info for buildvars file
env['buildvars'] = {}
# create list which will be populated with warnings if there are any
env['warnings'] = []

#################### Make sure install directories exist #####################

# This is for easybuild
if len(env['PREFIX']) != 0:
    env['prefix']=env['PREFIX']

if " " in env['prefix']:
    raise ValueError("Installation prefix contains space.")

env['BUILD_DIR'] = Dir(env['build_dir']).abspath
prefix = Dir(env['prefix']).abspath
env['buildvars']['prefix'] = prefix
env['incinstall'] = os.path.join(prefix, 'include')
env['bininstall'] = os.path.join(prefix, 'bin')
if IS_WINDOWS:
    env['libinstall'] = env['bininstall']
else:
    env['libinstall'] = os.path.join(prefix, 'lib','esys')
env['pyinstall']  = os.path.join(prefix, 'esys')
if env['trilinos_install'] == 'default':
    env['trilinos_install']=env['libinstall']

if not os.path.isdir(env['bininstall']):
    os.makedirs(env['bininstall'])
if not os.path.isdir(env['libinstall']):
    os.makedirs(env['libinstall'])
if not os.path.isdir(env['pyinstall']):
    os.makedirs(env['pyinstall'])

env.AppendUnique(CPPPATH = [env['incinstall']])
env.AppendUnique(LIBPATH = [env['libinstall']])


if not env['cxx_flags']:
    env['cxx_flags'] = env['cc_flags']
# default compiler/linker options
cxx_flags = [ '-std=c++20' ] # extra CXX flags (C++20 required for Trilinos 17/Kokkos 5)
cc_flags = [ ]
cc_optim = []
cc_debug = []
omp_flags = []
omp_ldflags = []
ld_flags = [ ]
fatalwarning = '' # switch to turn warnings into errors
sysheaderopt = env['sysheaderopt']  # how to indicate that a header is a system header

# env['CC'] might be a full path
cxx_name=os.path.basename(env['CXX'])
#
# this should be done via import of an appropriate template:

if cxx_name == 'icpc':
    # Intel compiler
    # #1478: class "std::auto_ptr<...>" was declared deprecated
    # #1875: offsetof applied to non-POD types is nonstandard (in boost)
    # removed -std=c99 because icpc doesn't like it and we aren't using c anymore
    cc_flags   += [ "-fPIC", "-w2", "-wd1875", "-wd1478", "-Wno-unknown-pragmas"]
    cc_optim    = [ "-Ofast", "-ftz", "-fno-alias", "-xCORE-AVX2", "-ipo"]
    #cc_optim    = "-Ofast -ftz -fno-alias -inline-level=2 -ipo -xCORE-AVX2"
    #cc_optim    = "-O2 -ftz -fno-alias -inline-level=2"
    #cc_optim    = "-O0 -ftz -fno-alias"
    #cc_optim    = "-O3 -ftz -fno-alias -inline-level=2 -ipo -xHost"
    cc_debug    = ["-g", "-O0", "-DDOASSERT", "-DDOPROF",  "-DBOUNDS_CHECK", "-DSLOWSHARECHECK"]
    omp_flags   = [ "-qopenmp" ]
    omp_ldflags = [ "-qopenmp" ] # removing -openmp-report (which is deprecated) because the replacement outputs to a file
    fatalwarning = [ "-Werror" ]
elif cxx_name.startswith('g++-mp'):
    # cc_flags +="-Wall" - switched off as there are problem in p4est compilation.
    cc_flags       += [ "-fPIC", "-fdiagnostics-color=always", "-Wno-uninitialized" ]
    cc_flags       += [ "-Wno-unknown-pragmas",  "-Wimplicit-function-declaration" ]
    cxx_flags      += ["-Wno-string-concatenation", "-Wno-unused-private-field" ]
    #if env['trilinos'] is True:
    #  cc_flags += "-Wno-unused-variable -Wno-exceptions -Wno-deprecated-declarations "
    cc_optim     = [ "-O3" ]
    cc_debug     = [ "-ggdb3","-O0", "-fdiagnostics-fixit-info", "-pedantic"]
    cc_debug    += ["-DDOASSERT" , "-DDOPROF", "-DBOUNDS_CHECK", "-DSLOWSHARECHECK" ]
    omp_flags    = ["-fopenmp"]
    omp_ldflags  = ["-fopenmp" ]
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"
elif cxx_name[:3] == 'g++':
    # GNU C++ on any system
    # note that -ffast-math is not used because it breaks isnan(),
    # see mantis #691
    cc_flags += ["-pedantic", "-Wall", "-fPIC", "-finline-functions" ]
    cc_flags += [ "-Wno-unknown-pragmas", "-Wno-sign-compare", "-Wno-system-headers", "-Wno-long-long", "-Wno-strict-aliasing" ]
    cc_flags += [ "-Wno-unused-function", "-Wno-narrowing" ]
    cc_flags += [ "-Wno-stringop-truncation", "-Wno-deprecated-declarations", "--param=max-vartrack-size=100000000"]
    cc_optim     = [ "-O3" ] # -march=native"
    #max-vartrack-size: avoid vartrack limit being exceeded with escriptcpp.cpp
    cc_debug     = [ "-g3", "-O0", "-DDOASSERT -DDOPROF", "-DBOUNDS_CHECK", "-DSLOWSHARECHECK", "--param=max-vartrack-size=100000000", '-D_GLIBCXX_DEBUG' ]
    omp_flags    = [ "-fopenmp"]
    omp_ldflags  = [ "-fopenmp"]
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"
elif cxx_name == 'cl':
    # Microsoft Visual C on Windows
    cc_flags     = ["/EHs", "/MD", "/GR", "/wd4068", "/D_USE_MATH_DEFINES", "/DDLL_HDF5"] # does this work
    cc_optim     = ["/O2", "/Op", "/W3"]
    cc_debug     = ["/Od", "/RTCcsu", "/ZI", "/DBOUNDS_CHECK"]
    fatalwarning = "/WX"
elif cxx_name == 'icl':
    # Intel C on Windows
    cc_flags     = ['/EHsc', '/GR', '/MD']
    cc_optim     = ['/fast', '/Oi', '/W3', '/Qssp', '/Qinline-factor-', '/Qinline-min-size=0', '/Qunroll']
    cc_debug     = ['/Od', '/RTCcsu', '/Zi', '/Y-', '/debug:all', '/Qtrapuv']
    omp_flags    = ['/Qvec-report0v/Qopenmp', '/Qopenmp-report0', '/Qparallel']
    omp_ldflags  = ['/Qvec-report0', '/Qopenmp', '/Qopenmp-report0', '/Qparallel']
elif cxx_name.startswith('clang++-mp'):
    # Clang++mp on any system
    # cc_flags +="-Wall" - switched off as there are problem in p4est compilation.
    cc_flags       += ["-fPIC", "-fdiagnostics-color=always", "-Wno-uninitialized", "-Wno-deprecated-declarations"]
    #cc_flags       += "-I/opt/local/lib/gcc12/gcc/arm64-apple-darwin22/12.2.0/include "
    cc_flags       += ["-Wno-unknown-pragmas" ]
    cxx_flags      += ["-Wno-string-concatenation", "-Wno-unused-private-field", "-Wimplicit-function-declaration" ]
    #if env['trilinos'] is True:
    #  cc_flags += "-Wno-unused-variable -Wno-exceptions -Wno-deprecated-declarations "
    cc_optim     = ["-O3" ]
    cc_debug     = ["-ggdb3 -O0", "-fdiagnostics-fixit-info", "-pedantic"]
    cc_debug    += ["-DDOASSERT", "-DDOPROF", "-DBOUNDS_CHECK", "-DSLOWSHARECHECK" ]
    omp_flags    = ["-fopenmp" ]
    omp_ldflags  = ["-fopenmp" ]
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"
elif cxx_name.startswith('clang++'):
    # Clang++ on any system (including Homebrew LLVM)
    cc_flags    += ["-Wall", "-fPIC", "-fdiagnostics-color=always",  "-Wno-uninitialized" ]
    cc_flags    += ["-Wno-unused-private-field", "-Wno-unknown-pragmas" ]
    cc_flags    += ["-Wno-deprecated-declarations"]
    #if env['trilinos'] is True:
    #  cc_flags += "-Wno-unused-variable -Wno-exceptions -Wno-deprecated-declarations "
    cc_optim     = ["-O3" ]  # -march=native"
    cc_debug     = ["-ggdb3", "-O0", "-fdiagnostics-fixit-info", "-pedantic" ]
    cc_debug    += ["-DDOASSERT", "-DDOPROF", "-DBOUNDS_CHECK", "-DSLOWSHARECHECK" ]
    omp_flags    = ["-fopenmp" ]
    omp_ldflags  = ["-fopenmp" ]
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"
elif cxx_name == 'mpic++':
    # MPIC++ wrapper - can wrap GCC or Clang depending on system
    # note that -ffast-math is not used because it breaks isnan()
    cc_flags += ["-pedantic", "-Wall", "-fPIC", "-finline-functions" ]
    cc_flags += ["-Wno-unknown-pragmas", "-Wno-sign-compare", "-Wno-system-headers", "-Wno-long-long", "-Wno-strict-aliasing" ]
    cc_flags += ["-Wno-unused-function", "-Wno-narrowing" ]
    cc_flags += ["-Wno-deprecated-declarations" ]
    # Removed GCC-specific flags: -Wno-stringop-truncation, --param=max-vartrack-size
    cc_optim     = ["-O3"]  # Removed -march=native which may not work on all systems
    cc_debug     = ["-g3", "-O0", "-DDOASSERT", "-DDOPROF", "-DBOUNDS_CHECK", "-DSLOWSHARECHECK" ]
    # Removed -D_GLIBCXX_DEBUG which is GCC-specific
    ld_flags += ["-fPIC", "-lmpi" ]
    # OpenMP library linking handled by -fopenmp linker flag (no explicit -lgomp/-lomp needed)
    omp_flags    = ["-fopenmp" ]
    omp_ldflags  = ["-fopenmp" ]
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"
elif cxx_name == 'mpiicpc':
    # note that -ffast-math is not used because it breaks isnan(),
    env['CC']= 'mpiicc'
    cxx_flags = ['-std=c++20' ]  # C++20 required for Trilinos 17/Kokkos 5
    cc_flags += ["-pedantic", "-Wall", "-fPIC", "-finline-functions"]
    cc_flags += ["-Wno-unknown-pragmas", "-Wno-sign-compare", "-Wno-system-headers", "-Wno-long-long", "-Wno-strict-aliasing"]
    cc_flags += ["-Wno-unused-function", "-Wno-narrowing"]
    cc_flags += ["-Wno-deprecated-declarations", "--param=max-vartrack-size=100000000"]
    cc_flags += ["-ipo" ]
    cc_optim     = ["-O2", "-march=native"]
    #max-vartrack-size: avoid vartrack limit being exceeded with escriptcpp.cpp
    cc_debug     = ["-g3", "-O0", "-DDOASSERT", "-DDOPROF", "-DBOUNDS_CHECK", "-DSLOWSHARECHECK", "--param=max-vartrack-size=100000000", '-D_GLIBCXX_DEBUG']
    ld_flags +=[ "-fPIC", "-lmpi" ]
    if env['openmp']:
      ld_flags += ["-lgomp" ]
    omp_flags    = ["-qopenmp" ]
    omp_ldflags  = ["-qopenmp" ]
    fatalwarning = "-Werror"
    sysheaderopt = "-isystem"


env['sysheaderopt']=sysheaderopt
# set defaults if not otherwise specified
if env['cc_flags']    == 'default': env['cc_flags'] = cc_flags
if env['cxx_flags']    == 'default': env['cxx_flags'] = cxx_flags
if env['cc_optim']    == 'default': env['cc_optim'] = cc_optim
if env['cc_debug']    == 'default': env['cc_debug'] = cc_debug
if env['omp_flags']   == 'default': env['omp_flags'] = omp_flags
if env['omp_ldflags'] == 'default': env['omp_ldflags'] = omp_ldflags
if env['ld_flags'] == 'default': env['ld_flags'] = ld_flags

if env['version_information'] != '0.0':
    env['SHLIBVERSION']=env['version_information']
    env['LDMODULEVERSION']=env['version_information']

if env['IS_OSX']:
    env.AppendUnique(SHLIBSONAMEFLAGS =  ["-Wl,-install_name=$_SHLIBSONAME"])
else:
    env.AppendUnique(SHLIBSONAMEFLAGS =   ["-Wl,-soname=$_SHLIBSONAME" ] )

if env['build_trilinos'] != 'never':
    if not os.path.isdir(env['trilinos_build']): # create a build folder if the user deleted it
        os.mkdir(env['trilinos_build'])

    # Sentinel file to track successful builds
    sentinel_file = os.path.join(env['trilinos_build'], '.trilinos_build_complete')

    # Determine if rebuild is needed
    needs_rebuild = False

    if env['build_trilinos'] == "always":
        needs_rebuild = True
        print("Trilinos: Always rebuild requested")
    elif not os.path.exists(sentinel_file):
        needs_rebuild = True
        print("Trilinos: No previous build found")
    elif env['build_trilinos'] in ["check", "make"]:
        # Check if source or config files are newer than the sentinel
        sentinel_time = os.path.getmtime(sentinel_file)

        # Check trilinos source directory - only check key CMake files for performance
        # If you modify source files, touch the top-level CMakeLists.txt to trigger rebuild
        trilinos_src_newer = False
        key_files = [
            os.path.join(env['trilinos_src'], 'CMakeLists.txt'),
            os.path.join(env['trilinos_src'], 'Version.cmake'),
        ]
        for fpath in key_files:
            if os.path.exists(fpath) and os.path.getmtime(fpath) > sentinel_time:
                trilinos_src_newer = True
                print(f"Trilinos: Key file {fpath} is newer than last build")
                break

        # Check build script files
        build_scripts = []
        if env['trilinos_make_sh'] == 'default':
            if env['mpi'] != 'none':
                build_scripts.append('scripts/trilinos_mpi.sh')
            else:
                build_scripts.append('scripts/trilinos_nompi.sh')
        else:
            build_scripts.append(env['trilinos_make_sh'])

        scripts_newer = any(os.path.getmtime(s) > sentinel_time for s in build_scripts if os.path.exists(s))

        # Check if esys.trilinos directory is missing
        trilinos_install_missing = not os.path.isdir(env['trilinos_install'])

        needs_rebuild = trilinos_src_newer or scripts_newer or trilinos_install_missing

        if trilinos_src_newer:
            print("Trilinos: Source files changed")
        if scripts_newer:
            print("Trilinos: Build scripts changed")
        if trilinos_install_missing:
            print("Trilinos: Installation directory missing")
        if not needs_rebuild:
            print("Trilinos: Up to date, skipping build")

    if needs_rebuild:
        if env['openmp']:
            OPENMPFLAG='ON'
        else:
            OPENMPFLAG='OFF'
        if not env['cc'] == 'default ':
            os.environ['CC'] = env['cc']
        if not env['cxx'] == 'default ':
            os.environ['CXX'] = env['cxx']
        # Prepare METIS and ParMETIS info for Trilinos build
        # Extract paths from prefix (can be single path or [inc_path, lib_path])
        def get_inc_lib_paths(prefix):
            if isinstance(prefix, (list, tuple)) and len(prefix) >= 2:
                return prefix[0], prefix[1]
            elif isinstance(prefix, (list, tuple)) and len(prefix) == 1:
                return prefix[0], prefix[0]
            else:
                return str(prefix), str(prefix)

        metis_enable = 'ON' if env['metis'] else 'OFF'
        metis_inc, metis_lib = get_inc_lib_paths(env['metis_prefix']) if env['metis'] else ('', '')
        parmetis_enable = 'ON' if env['parmetis'] else 'OFF'
        parmetis_inc, parmetis_lib = get_inc_lib_paths(env['parmetis_prefix']) if env['parmetis'] else ('', '')
        scotch_enable = 'ON' if env['scotch'] else 'OFF'
        scotch_inc, scotch_lib = get_inc_lib_paths(env['scotch_prefix']) if env['scotch'] else ('', '')

        SHARGS = [ env['trilinos_install'], env['CC'],  env['CXX'], OPENMPFLAG, env['trilinos_src'],
                   metis_enable, metis_inc, metis_lib, parmetis_enable, parmetis_inc, parmetis_lib,
                   scotch_enable, scotch_inc, scotch_lib ]

        print("Initialization of Trilinos build using", SHARGS )
        if env['trilinos_make_sh'] == 'default':
            if env['mpi'] != 'none':
                shutil.copy("scripts/trilinos_mpi.sh", env['trilinos_build'] + "/trilinos_mpi.sh")
                print("Building (MPI) trilinos..............................")
                SH ="trilinos_mpi.sh"
            else:
                shutil.copy("scripts/trilinos_nompi.sh", env['trilinos_build'] + "/trilinos_nompi.sh")
                print("Building (no MPI) trilinos..............................")
                SH = "trilinos_nompi.sh"
        else:
            shutil.copy(env['trilinos_make_sh'], os.path.join(env['trilinos_build'],"hostmake.sh"))
            SH = "hostmake.sh"
        import subprocess
        p_init = subprocess.run(['sh', SH ] + SHARGS, capture_output=False, cwd =  env['trilinos_build'],  text=True)
        if p_init.returncode :
            print(">>> Initialization of Trilinos build failed. Scons stopped.")
            Exit(1)

        SHARGS = [ ]
        if env['build_trilinos'] == "always" :
                SHARGS += ['--always-make']
        SHARGS += [ f'-j{GetOption("num_jobs")}', 'install' ]
        print("Trilinos build using", SHARGS)
        p_make = subprocess.run( [ 'make' ] + SHARGS, capture_output=False, cwd =  env['trilinos_build'],  text=True)
        if p_make.returncode :
            print(">>> Installation of Trilinos failed. Scons stopped.")
            Exit(1)

        # Create/update sentinel file on successful build
        with open(sentinel_file, 'w') as f:
            import time
            f.write(f"Trilinos build completed at {time.ctime()}\n")
        print(f"Trilinos: Build complete, sentinel file updated")

    env['trilinos'] = True
    env['trilinos_prefix'] = env['trilinos_install']
#=====
if env['longindices']:
  if env['paso']:
    env.Append(CPPDEFINES = ['ESYS_INDEXTYPE_LONG'])
  else:
    env['warnings'].append("The longindices feature requires paso!")

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

# --- set  options:
if env['cc_flags'] != 'default': env.Append(CCFLAGS = env['cc_flags'])
if env['cxx_flags'] != 'default': env.Append(CXXFLAGS = env['cxx_flags'])
if env['debug']:
    if env['cc_debug'] != 'default':
        env.AppendUnique(CCFLAGS = env['cc_debug'])
        env.AppendUnique(CXXFLAGS = env['cc_debug'])
else:
    if env['cc_optim'] != 'default':
        env.AppendUnique(CCFLAGS = env['cc_optim'])
        env.AppendUnique(CXXFLAGS = env['cc_optim'])
env['buildvars']['debug']=int(env['debug'])
if env['ld_flags'] != 'default': env.Append(LINKFLAGS =env['ld_flags'])

if env['cxx_extra'] : env.AppendUnique(CXXFLAGS = env['cxx_extra'])
if env['cc_extra'] : env.AppendUnique(CCFLAGS = env['cc_extra'])
if env['ld_extra']  : env.AppendUnique(LINKFLAGS =env['ld_extra'])


# Disable OpenMP if no flags provided
if env['openmp'] and env['omp_flags'] == '':
   env['warnings'].append("OpenMP requested but no flags provided - disabling OpenMP!")
   env['openmp'] = False

if env['openmp']:
    env.AppendUnique(CCFLAGS = env['omp_flags'])
    if env['omp_ldflags']: env.AppendUnique(LINKFLAGS = env['omp_ldflags'])
else:
    env['omp_flags']= [ ]
    env['omp_ldflags']= [ ]

# Flags used by the p4est program
if env['openmp']:
    env.Append(CPPDEFINES=['P4EST_HAVE_OPENMP'])
    env.Append(CPPDEFINES=['P4EST_OPENMP'])
    env.Append(CPPDEFINES=['P4EST_ENABLE_OPENMP'])
    env.Append(CPPDEFINES=['ENABLE_OPENMP'])
    env.Append(CPPDEFINES=['HAVE_OPENMP'])
    env.Append(CPPDEFINES=['OPENMP'])

env['buildvars']['openmp']=int(env['openmp'])
env.AppendUnique(LIBS = env['sys_libs'])


global_revision=''
try:
    global_revision = os.popen('git show -s --format=%ct').read()
    global_revision = re.sub(':.*', '', global_revision)
    global_revision = re.sub('[^0-9]', '', global_revision)
    print("Got git revision ", global_revision)
except:
    env['warnings'].append("Could not detect the git commit information!")
env['svn_revision']=global_revision
env['buildvars']['svn_revision']=global_revision
env.Append(CPPDEFINES=['GIT_BUILD'])



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
    sl_suffix = '.pyd' if IS_WINDOWS else '.so'
    return env.SharedLibrary(target, source, SHLIBPREFIX='', SHLIBSUFFIX=sl_suffix)
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

######## optional python modules (sympy, pyproj)
env=checkOptionalModules(env)

######## optional dependencies (HDF5, MKL, UMFPACK, MUMPS, Lapack, Silo, ...)

env=checkOptionalLibraries(env)

######## PDFLaTeX (for documentation)
env=checkPDFLatex(env)


# =================================
# set defaults for launchers if not otherwise specified
if env['prelaunch'] == 'default':
    if env['mpi'] == 'INTELMPI' and env['openmp']:
        env['prelaunch'] = "export I_MPI_PIN_DOMAIN=omp"
    # elif env['mpi'] == 'OPENMPI':
        # transform comma-separated list to '-x a -x b -x c ...'
        # env['prelaunch'] = "EE=$echo -x %e|sed -e 's/,/ -x /g'"
    elif env['mpi'] == 'MPT':
        env['prelaunch'] = "export MPI_NUM_MEMORY_REGIONS=0"
    elif env['mpi'] == 'MPICH2':
        env['prelaunch'] = "mpdboot -n %n -r ssh -f %f"
    else:
        env['prelaunch'] = ""

# Used by p4est
if env['mpi'] != 'none':
    env.Append(CPPDEFINES = ['P4EST_ENABLE_MPI'])
    env.Append(CPPDEFINES = ['P4EST_ENABLE_MPICOMMSHARED'])
    env.Append(CPPDEFINES = ['P4EST_ENABLE_MPIIO'])
    env.Append(CPPDEFINES = ['P4EST_ENABLE_MPISHARED'])
    env.Append(CPPDEFINES = ['P4EST_ENABLE_MPITHREAD'])
    env.Append(CPPDEFINES = ['P4EST_ENABLE_MPIWINSHARED'])
    env.Append(CPPDEFINES = ['P4EST_MPI'])
    env.Append(CPPDEFINES = ['P4EST_MPIIO'])
    env.Append(CPPDEFINES = ['SC_ENABLE_MPI'])
    env.Append(CPPDEFINES = ['SC_ENABLE_MPICOMMSHARED'])
    env.Append(CPPDEFINES = ['SC_ENABLE_MPIIO'])
    env.Append(CPPDEFINES = ['SC_ENABLE_MPISHARED'])
    env.Append(CPPDEFINES = ['SC_ENABLE_MPITHREAD'])
    env.Append(CPPDEFINES = ['SC_ENABLE_MPIWINSHARED'])
    env.Append(CPPDEFINES = ['SC_MPI'])
    env.Append(CPPDEFINES = ['SC_MPIIO'])

if env['launcher'] == 'default':
    if env['mpi'] == 'INTELMPI':
        env['launcher'] = "mpirun -hostfile %f -n %N -ppn %p %b"
    elif env['mpi'] == 'OPENMPI':
        if env['mpi_no_host']:
            hostoptionstr=''
        else:
            hostoptionstr='--host %h'
        # default to OpenMPI version 1.10 or higher
        env['launcher'] = "mpirun ${AGENTOVERRIDE} --gmca mpi_warn_on_fork 0 ${EE} "+hostoptionstr+" --map-by node:pe=%t -bind-to core -np %N %b"
        if 'orte_version' in env:
            major,minor,point = [int(i) for i in env['orte_version'].split('.')]
            if major == 1 and minor < 10:
                env['launcher'] = "mpirun ${AGENTOVERRIDE} --gmca mpi_warn_on_fork 0 ${EE} "+hostoptionstr+" --cpus-per-rank %t -np %N %b"
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
else:
   print("Building escript with the domains:  %s"%(", ".join(env['domains'])))

# keep some of our install paths first in the list for the unit tests
env.PrependENVPath(LD_LIBRARY_PATH_KEY, env['libinstall'])
env.PrependENVPath('PYTHONPATH', prefix)
env['ENV']['ESCRIPT_ROOT'] = prefix

if not env['verbose']:
    env['CCCOMSTR'] = "Compiling $TARGET"
    env['SHCCCOMSTR'] = "Compiling $TARGET"
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
    env.AppendUnique(CCFLAGS = fatalwarning)

Export(
  ['env',
   'dodgy_env',
   'IS_WINDOWS', 'IS_OSX',
   'TestGroups'
  ]
)

target_init = env.Command(os.path.join(env['pyinstall'],'__init__.py'), None, Touch('$TARGET'))
env.Alias('target_init', [target_init])

# escript can't be turned off
build_all_list = ['build_escript']
install_all_list = ['target_init', 'install_escript']

#p4est
if "oxley" in env['domains'] and not env['zlib']:
    print("Error: oxley domain requires zlib. Please set zlib=True in your options file.")
    print("On Debian/Ubuntu, install with: sudo apt-get install zlib1g-dev")
    env.Exit(1)

if env['p4est'] and "oxley" in env['domains']  :
    build_all_list += ['build_p4est']
    install_all_list += ['install_p4est']
    env['p4est']=True
    env['p4est_libs']=['p4est','sc']
    env['escript_src']=os.getcwd()

if env['usempi']:
    build_all_list += ['build_pythonMPI', 'build_overlord']
    install_all_list += ['install_pythonMPI', 'install_overlord']

env['buildvars']['paso'] = int(env['paso'])
if env['paso']:
    env.Append(CPPDEFINES = ['ESYS_HAVE_PASO'])
    build_all_list += ['build_paso']
    install_all_list += ['install_paso']

env['buildvars']['trilinos'] = int(env['trilinos'])
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
    if 'finley' in env['domains']:
        build_all_list += ['build_escriptreader']
        install_all_list += ['install_escriptreader']




variant='$BUILD_DIR/$PLATFORM/'
env.SConscript('escriptcore/SConscript', variant_dir=variant+'escriptcore', duplicate=0)
env.SConscript('escript/SConscript', variant_dir=variant+'escript', duplicate=0)
env.SConscript('pythonMPI/SConscript', variant_dir=variant+'pythonMPI', duplicate=0)
env.SConscript('tools/overlord/SConscript', variant_dir=variant+'tools/overlord', duplicate=0)
env.SConscript('paso/SConscript', variant_dir=variant+'paso', duplicate=0)
env.SConscript('trilinoswrap/SConscript', variant_dir=variant+'trilinoswrap', duplicate=0)
env.SConscript('cusplibrary/SConscript')
env.SConscript('finley/SConscript', variant_dir=variant+'finley', duplicate=0)
if env['p4est']:
    env.SConscript('p4est/SConscript', variant_dir=variant+'p4est', duplicate=0)
if os.path.isdir('oxley'):
    env.SConscript('oxley/SConscript', variant_dir=variant+'oxley', duplicate=0)
env.SConscript('ripley/SConscript', variant_dir=variant+'ripley', duplicate=0)
env.SConscript('speckley/SConscript', variant_dir=variant+'speckley', duplicate=0)
env.SConscript('weipa/SConscript', variant_dir=variant+'weipa', duplicate=0)
env.SConscript('tools/escriptconvert/SConscript', variant_dir=variant+'tools/escriptconvert', duplicate=0)
env.SConscript('doc/SConscript', variant_dir=variant+'doc', duplicate=0)

env.Alias('build', build_all_list)

install_all_list += [env.Install(Dir('scripts',env['build_dir']), os.path.join('scripts', 'release_sanity.py'))]

if env['mpi']:
    install_all_list += ['install_pythonMPI']

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
if env['domains'] == all_domains and not env['insane']:
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

env.Alias('pdfdocs',['user_pdf'  ])
env.Alias('basedocs', ['pdfdocs','examples_tarfile', 'examples_zipfile', 'api_doxygen', 'readme_html', 'installation_html'])
env.Alias('docs', ['basedocs', 'sphinxdoc'])
env.Alias('release_prep', ['docs', 'install'])
# Epydoc removed - use Sphinx instead
#env.Alias('release_prep_old', ['basedocs', 'api_epydoc', 'install'])

# The test scripts are always generated, this target allows us to
# generate the testscripts without doing a full build
env.Alias('testscripts',[])

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
    print(f"*** Config Summary (see config.log and {env['libinstall']}/buildvars for details) ***")
    print("Escript revision %s"%global_revision)
    print("  Install prefix:  %s"%env['prefix'])
    print("          Python:  %s (Version %s)"%(env['pythoncmd'],env['python_version']))
    print("           boost:  %s (Version %s)"%(env['boost_prefix'],env['boost_version']))
    if env['have_boost_numpy']:
        print("     boost numpy:  YES")
    else:
        print("     boost numpy:  NO")
    if env['build_trilinos']:
        print("        trilinos:  %s (built-in, Version %s)" % (env['trilinos_prefix'], env['trilinos_version']))
    else:    
        if env['trilinos']:
            print("        trilinos:  %s (Version %s)" % (env['trilinos_prefix'],env['trilinos_version']))
        else:
            print("        trilinos:  NO")
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
    if env['mpi4py']:
        print("          mpi4py:  YES")
    else:
        print("          mpi4py:  NO")
    if env['metis']:
        print("           METIS:  %s (Version %s)"%(env['metis_prefix'],env['metis_version']))
    else:
        d_list.append('metis')
    if env['parmetis']:
        print("        ParMETIS:  %s (Version %s)"%(env['parmetis_prefix'],env['parmetis_version']))
    else:
        d_list.append('parmetis')
    if env['scotch']:
        print("          Scotch:  %s (Version %s)"%(env['scotch_prefix'],env['scotch_version']))
    else:
        d_list.append('scotch')
    if env['uselapack']:
        print("          LAPACK:  YES (flavour: %s)"%env['lapack'])
    else:
        d_list.append('lapack')
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
        if env['mumps_seq']:
            direct.append('mumps')
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
    for i in ('hdf5', 'weipa','debug','openmp','cppunit','mkl','mpi4py', 'zlib',
             'mumps_seq', 'netcdf', 'scipy','silo','sympy','umfpack','visit'):
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

    print("Add to PYTHONPATH :", env['prefix'])
    print("Add to LD_LIBRARY_PATH :", env['libinstall'],":",os.path.join(env['trilinos_install'],'lib'))
    print("If publishing results using esys-escript, please cite us: see https://github.com/LutzGross/esys-escript.github.io/")

atexit.register(print_summary)
