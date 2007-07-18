#         Copyright 2006 by ACcESS MNRF
#
#              http://www.access.edu.au
#       Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php

# top-level Scons configuration file for all esys13 modules
# Begin initialisation Section
# all of this section just intialises default environments and helper
# scripts. You shouldn't need to modify this section.
EnsureSConsVersion(0,96,91)
EnsurePythonVersion(2,3)

#===============================================================
#   import tools:
import glob
import sys, os
import socket
# Add our extensions
if sys.path.count('scons')==0: sys.path.append('scons')
import scons_extensions

#===============================================================
#   check on windows or linux platform
#
IS_WINDOWS_PLATFORM = (os.name== "nt")

if IS_WINDOWS_PLATFORM:
   tools_prefix="C:\\Program Files\\"
else:
   tools_prefix="/usr"

#==============================================================================================     
#    
#    get the installation prefix
#
prefix = ARGUMENTS.get('prefix', Dir('#.').abspath)

# We may also need to know where python's site-packages subdirectory lives
python_version = 'python%s.%s'%(sys.version_info[0],sys.version_info[1])

if prefix == "/usr":
   # Install as a standard python package in /usr/lib64 if available, else in /usr/lib
   if os.path.isdir(  prefix+"/lib64/"+python_version+"/site-packages"):
      dir_packages =  prefix+"/lib64/"+python_version+"/site-packages"
      dir_libraries = prefix+"/lib64"
   elif os.path.isdir(prefix+"/lib/"+python_version+"/site-packages"):
      dir_packages =  prefix+"/lib/"+python_version+"/site-packages"
      dir_libraries = prefix+"/lib"
   else:
      print "Install prefix is /usr but couldn't find python package directory in either"
      print "/usr/lib64/"+python_version+"/site-packages or /usr/lib/"+python_version+"/site-packages"
      sys.exit(1)
   dir_examples = prefix+"/share/doc/esys"
else:
   # Install using the usual escript directory structure
   dir_packages = prefix
   dir_libraries = prefix+"/lib"
   dir_examples = prefix
dir_packages += "/esys"
dir_examples += "/examples"

print "Install prefix is: ", prefix
print "	python packages will be installed in:	", dir_packages
print "	libraries will be installed in:		", dir_libraries
print "	examples will be installed in:		", dir_examples

#==============================================================================================     

# Default options and options help text
# These are defaults and can be overridden using command line arguments or an options file.
# if the options_file or ARGUMENTS do not exist then the ones listed as default here are used
# DO NOT CHANGE THEM HERE
# Where to install?
#==============================================================================================     
#    
#    get the options file if present:
#
if ARGUMENTS.get('options_file',0):
   options_file = ARGUMENTS.get('options_file',0)
else:
   from string import ascii_letters,digits
   hostname=""
   for s in socket.gethostname().split('.')[0]:
      if s in ascii_letters+digits:
         hostname+=s
      else:
         hostname+="_"
   options_file = os.path.join("scons",hostname+"_options.py")
   
if os.path.isfile(options_file):
   print "option file is ",options_file,"."
else:
   print "option file is ",options_file, "(not present)."
# and load it
opts = Options(options_file, ARGUMENTS)
#================================================================
#
#   check if UMFPACK is installed on the system:
#
uf_root=None
for i in [ 'UMFPACK', 'umfpack', 'ufsparse', 'UFSPARSE']:
   if os.path.isdir(os.path.join(tools_prefix,'include',i)): 
       uf_root=i
       print i," is used form ",tools_prefix
       break
if not uf_root==None:
   umf_path_default=os.path.join(tools_prefix,'include',uf_root)
   umf_lib_path_default=os.path.join(tools_prefix,'lib')
   umf_libs_default=['umfpack']
   amd_path_default=os.path.join(tools_prefix,'include',uf_root)
   amd_lib_path_default=os.path.join(tools_prefix,'lib')
   amd_libs_default=['amd']
   ufc_path_default=os.path.join(tools_prefix,'include',uf_root)
else:
   umf_path_default=None
   umf_lib_path_default=None
   umf_libs_default=None
   amd_path_default=None
   amd_lib_path_default=None
   amd_libs_default=None
   ufc_path_default=None
#
#==========================================================================
#
#    python installation:
#
if IS_WINDOWS_PLATFORM:
   python_path_default=os.path.join(tools_prefix,'python%s%s'%(sys.version_info[0],sys.version_info[1]),"include")
   python_lib_path_default=os.path.join(tools_prefix,'python%s%s'%(sys.version_info[0],sys.version_info[1]),"libs")
   python_libs_default=["python%s%s"%(sys.version_info[0],sys.version_info[1])]
else:
   python_path_default=os.path.join(tools_prefix,'include','python%s.%s'%(sys.version_info[0],sys.version_info[1]))
   python_lib_path_default=os.path.join(tools_prefix,'lib')
   python_libs_default=["python%s.%s"%(sys.version_info[0],sys.version_info[1])]

#==========================================================================
#
#    boost installation:
#
if IS_WINDOWS_PLATFORM:
   boost_libs_path_default=os.path.join(tools_prefix,'boost','lib')
   boost_libs_default=None
   for i in os.listdir(boost_libs_path_default): 
      name=os.path.splitext(i)
      if name[1] == ".dll" and name[0].startswith("boost_python"):
          if boost_libs_default == None: 
	     boost_libs_default= [ name[0] ] 
	  else:
	     if not name[0].find("-gd-"): boost_libs_default=[ name[0] ]
   boost_path_default=os.path.join(tools_prefix,'boost','include','boost-%s'%(boost_libs_default[0].split("-")[-1],))
else:
   boost_path_default=os.path.join(tools_prefix,'include')
   boost_libs_path_default=os.path.join(tools_prefix,'lib')
   boost_libs_default=['boost_python']
#==========================================================================
#
#    check if netCDF is installed on the system:
#
if IS_WINDOWS_PLATFORM:
    netcdf_dir=os.path.join(tools_prefix,'netcdf')
    netCDF_path_default=os.path.join(netcdf_dir,'include')
    netCDF_lib_path_default=os.path.join(netcdf_dir,'lib')
else:
    netCDF_path_default=os.path.join(tools_prefix,'include','netcdf-3')
    netCDF_lib_path_default=os.path.join(tools_prefix,'lib','netcdf-3')

if os.path.isdir(netCDF_path_default) and os.path.isdir(netCDF_lib_path_default):
     useNetCDF_default='yes'
     netCDF_libs_default=[ 'netcdf_c++', 'netcdf' ]
else:
     useNetCDF_default='no'
     netCDF_path_default=None
     netCDF_lib_path_default=None
     netCDF_libs_default=None

if IS_WINDOWS_PLATFORM: 
        useNetCDF_default='no' # be default netcdf is not supported on windows. 
#==========================================================================
#
#    compile:
#
if IS_WINDOWS_PLATFORM:
    # cc_flags_default  = '/GR /EHsc /MD /Qc99 /Qopenmp /Qopenmp-report1 /O3 /G7 /Qprec /Qpar-report1 /QxP /QaxP'
    # cc_flags_debug_default  = '/Od /MDd /RTC1 /GR /EHsc /Qc99 /Qopenmp /Qopenmp-report1 /Qprec'
    cc_flags_default  = '/nologo /EHsc /GR  /wd4068 /O2 /Op /MT /W3 /Ob0 /Z7'
    cc_flags_debug_default  ='/nologo /EHsc /GR  /wd4068 /Od /RTC1 /MTd /ZI /Ob0 /Z7'
    
    cc_flags_default  = '/nologo /EHsc /GR  /O2 /MT /W3 /Ob0 /Z7 /wd4068'
    cc_flags_debug_default  ='/nologo /EHsc /GR /Od /RTC1 /MTd /W3 /Ob0 /Z7/wd4068'
    cxx_flags_default = ''
    cxx_flags_debug_default = ''
    cc_common_flags = '/FD /EHsc /GR /wd4068 '
else:
   cc_flags_default='-O3 -std=c99 -ffast-math -fpic -Wno-unknown-pragmas -ansi -pedantic-errors'
   cc_flags_debug_default='-g -O0 -ffast-math -std=c99 -fpic -Wno-unknown-pragmas -ansi -pedantic-errors'
   cxx_flags_default='--no-warn -ansi'
   cxx_flags_debug_default='--no-warn -ansi -DDOASSERT'
#==============================================================================================     
# Default options and options help text
# These are defaults and can be overridden using command line arguments or an options file.
# if the options_file or ARGUMENTS do not exist then the ones listed as default here are used
# DO NOT CHANGE THEM HERE
opts.AddOptions(
# Where to install esys stuff
  ('incinstall', 'where the esys headers will be installed',             Dir('#.').abspath+'/include'),
  ('libinstall', 'where the esys libraries will be installed',           dir_libraries),
  ('pyinstall', 'where the esys python modules will be installed',       dir_packages),
  ('exinstall', 'where the esys examples will be installed',             dir_examples),
  ('src_zipfile', 'the source zip file will be installed.',              Dir('#.').abspath+"/release/escript_src.zip"),
  ('test_zipfile', 'the test zip file will be installed.',               Dir('#.').abspath+"/release/escript_tests.zip"),
  ('src_tarfile', 'the source tar file will be installed.',              Dir('#.').abspath+"/release/escript_src.tar.gz"),
  ('test_tarfile', 'the test tar file will be installed.',               Dir('#.').abspath+"/release/escript_tests.tar.gz"),
  ('examples_tarfile', 'the examples tar file will be installed.',       Dir('#.').abspath+"/release/doc/escript_examples.tar.gz"),
  ('examples_zipfile', 'the examples zip file will be installed.',       Dir('#.').abspath+"/release/doc/escript_examples.zip"),
  ('guide_pdf', 'name of the user guide in pdf format',                  Dir('#.').abspath+"/release/doc/user/guide.pdf"),
  ('api_epydoc', 'name of the epydoc api docs directory',                Dir('#.').abspath+"/release/doc/epydoc"),
  ('guide_html', 'name of the directory for user guide in html format',  Dir('#.').abspath+"/release/doc/user/html"),
  ('api_doxygen', 'name of the doxygen api docs directory',prefix+"/release/doc/doxygen"),
# Compilation options
  BoolOption('dodebug', 'Do you want a debug build?', 'no'),
  ('options_file', "Optional file containing preferred options. Ignored if it doesn't exist (default: scons/<hostname>_options.py)", options_file),
  ('cc_defines','C/C++ defines to use', None),
  ('cc_flags','C compiler flags to use (Release build)', cc_flags_default),
  ('cc_flags_debug', 'C compiler flags to use (Debug build)', cc_flags_debug_default),
  ('cxx_flags', 'C++ compiler flags to use (Release build)', cxx_flags_default),
  ('cxx_flags_debug', 'C++ compiler flags to use (Debug build)', cxx_flags_debug_default),
  ('ar_flags', 'Static library archiver flags to use', None),
  ('sys_libs', 'System libraries to link with', None),
  ('tar_flags','flags for zip files','-c -z'),
# MKL
  PathOption('mkl_path', 'Path to MKL includes', None),
  PathOption('mkl_lib_path', 'Path to MKL libs', None),
  ('mkl_libs', 'MKL libraries to link with', None),
# SCSL
  PathOption('scsl_path', 'Path to SCSL includes', None),
  PathOption('scsl_lib_path', 'Path to SCSL libs', None),
  ('scsl_libs', 'SCSL libraries to link with', None),
  ('scsl_libs_MPI', 'SCSL libraries to link with for MPI build', None),
# UMFPACK
  PathOption('ufc_path', 'Path to UFconfig includes', ufc_path_default),
  PathOption('umf_path', 'Path to UMFPACK includes', umf_path_default),
  PathOption('umf_lib_path', 'Path to UMFPACK libs', umf_lib_path_default),
  ('umf_libs', 'UMFPACK libraries to link with', umf_libs_default),
# AMD (used by UMFPACK)
  PathOption('amd_path', 'Path to AMD includes', amd_path_default),
  PathOption('amd_lib_path', 'Path to AMD libs', amd_lib_path_default),
  ('amd_libs', 'AMD libraries to link with', amd_libs_default),
# BLAS
  PathOption('blas_path', 'Path to BLAS includes', None),
  PathOption('blas_lib_path', 'Path to BLAS libs', None),
  ('blas_libs', 'BLAS libraries to link with', None),
# netCDF
  ('useNetCDF', 'switch on/off the usage of netCDF', useNetCDF_default),
  PathOption('netCDF_path', 'Path to netCDF includes', netCDF_path_default),
  PathOption('netCDF_lib_path', 'Path to netCDF libs', netCDF_lib_path_default),
  ('netCDF_libs', 'netCDF C++ libraries to link with', netCDF_libs_default),
# Python
# locations of include files for python
  PathOption('python_path', 'Path to Python includes', python_path_default),
  PathOption('python_lib_path', 'Path to Python libs', python_lib_path_default),
  ('python_libs', 'Python libraries to link with', python_libs_default),
# Boost
  PathOption('boost_path', 'Path to Boost includes', boost_path_default),
  PathOption('boost_libs_path', 'Path to Boost libs', boost_libs_path_default),
  ('boost_libs', 'Boost libraries to link with', boost_libs_default),
# Doc building
#  PathOption('doxygen_path', 'Path to Doxygen executable', None),
#  PathOption('epydoc_path', 'Path to Epydoc executable', None),
# PAPI
  PathOption('papi_path', 'Path to PAPI includes', None),
  PathOption('papi_lib_path', 'Path to PAPI libs', None),
  ('papi_libs', 'PAPI libraries to link with', None),
# MPI
  BoolOption('useMPI', 'Compile parallel version using MPI', 'no'),
)
#=================================================================================================
#
#   Note: On the Altix the intel compilers are not automatically
#   detected by scons intelc.py script. The Altix has a different directory
#   path and in some locations the "modules" facility is used to support
#   multiple compiler versions. This forces the need to import the users PATH
#   environment which isn't the "scons way"
#   This doesn't impact linux and windows which will use the default compiler (g++ or msvc, or the intel compiler if it is installed on both platforms)
#   FIXME: Perhaps a modification to intelc.py will allow better support for ia64 on altix
#
if IS_WINDOWS_PLATFORM:
      env = Environment(tools = ['default', 'msvc'], options = opts)
else:
   if os.uname()[4]=='ia64':
      env = Environment(tools = ['default', 'intelc'], options = opts)
      if env['CXX'] == 'icpc':
         env['LINK'] = env['CXX'] # version >=9 of intel c++ compiler requires use of icpc to link in C++ runtimes (icc does not). FIXME: this behaviour could be directly incorporated into scons intelc.py
   else:
      env = Environment(tools = ['default'], options = opts)
Help(opts.GenerateHelpText(env))
#=================================================================================================
#
#     Initialise Scons Build Environment
#     check for user environment variables we are interested in
try:
   python_path = os.environ['PYTHONPATH']
   env['ENV']['PYTHONPATH'] = python_path
except KeyError:
   python_path = ''

try:
   omp_num_threads = os.environ['OMP_NUM_THREADS']
except KeyError:
   omp_num_threads = 1
env['ENV']['OMP_NUM_THREADS'] = omp_num_threads

try:
   env['ENV']['DISPLAY'] = os.environ['DISPLAY']
   env['ENV']['XAUTHORITY'] = os.environ['XAUTHORITY']
except KeyError:
   pass

try:
   path = os.environ['PATH']
   env['ENV']['PATH'] = path
except KeyError:
   path = ''
try:
   ld_library_path = os.environ['LD_LIBRARY_PATH']
   env['ENV']['LD_LIBRARY_PATH'] = ld_library_path
except KeyError:
   ld_library_path = ''
#==========================================================================
#
#    Add some customer builders
#
py_builder = Builder(action = scons_extensions.build_py, suffix = '.pyc', src_suffix = '.py', single_source=True)
env.Append(BUILDERS = {'PyCompile' : py_builder});

if IS_WINDOWS_PLATFORM:
   runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed', src_suffix='.exe', single_source=True)
else:
   runUnitTest_builder = Builder(action = scons_extensions.runUnitTest, suffix = '.passed', single_source=True)
env.Append(BUILDERS = {'RunUnitTest' : runUnitTest_builder});

runPyUnitTest_builder = Builder(action = scons_extensions.runPyUnitTest, suffix = '.passed', src_suffic='.py', single_source=True)
env.Append(BUILDERS = {'RunPyUnitTest' : runPyUnitTest_builder});

# Convert the options which are held in environment variable into python variables for ease of handling and configure compilation options
try:
   incinstall = env['incinstall']
   env.Append(CPPPATH = [incinstall,])
except KeyError:
   incinstall = None
try:
   libinstall = env['libinstall']
   env.Append(LIBPATH = [libinstall,]) # Adds -L for building of libescript.so libfinley.so escriptcpp.so finleycpp.so
   env.PrependENVPath('LD_LIBRARY_PATH', libinstall)
   if env['PLATFORM'] == "win32":
      env.PrependENVPath('PATH', libinstall)
      env.PrependENVPath('PATH', env['boost_libs_path'])
except KeyError:
   libinstall = None
try:
   pyinstall = env['pyinstall'] # all targets will install into pyinstall/esys but PYTHONPATH points at straight pyinstall so you go import esys.escript etc
   env.PrependENVPath('PYTHONPATH', env['pyinstall'])
except KeyError:
   pyinstall = None
try:
   exinstall = env['exinstall']
except KeyError:
   exinstall = None
try:
   dodebug = env['dodebug']
except KeyError:
   dodebug = None
try:
   useMPI = env['useMPI']
except KeyError:
   useMPI = None
try:
   cc_defines = env['cc_defines']
   env.Append(CPPDEFINES = cc_defines)
except KeyError:
   pass

if dodebug:
  if useMPI:
    try:
      flags = env['cc_flags_debug_MPI']
      env.Append(CCFLAGS = flags)
    except KeyError:
      pass
  else:
    try:
      flags = env['cc_flags_debug']
      env.Append(CCFLAGS = flags)
    except KeyError:
      pass
else:
  if useMPI:
   try:
     flags = env['cc_flags_MPI']
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
   if useMPI:
     try:
        flags = env['cxx_flags_debug_MPI']
        env.Append(CXXFLAGS = flags)
     except KeyError:
        pass
   else:
     try:
        flags = env['cxx_flags_debug']
        env.Append(CXXFLAGS = flags)
     except KeyError:
        pass
else:
   if useMPI:
     try:
        flags = env['cxx_flags_MPI']
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
   sys_libs = []

try:
   tar_flags = env['tar_flags']
   env.Replace(TARFLAGS = tar_flags)
except KeyError:
   pass

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

if useMPI:
   mkl_libs = []
else:
   try:
      mkl_libs = env['mkl_libs']
   except KeyError:
      mkl_libs = []

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

if useMPI:
  try:
    scsl_libs = env['scsl_libs_MPI']
  except KeyError:
    scsl_libs = []
else:
  try:
    scsl_libs = env['scsl_libs']
  except KeyError:
    scsl_libs = []

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

if useMPI:
    umf_libs = []
else:
   try:
      umf_libs = env['umf_libs']
   except KeyError:
      umf_libs = []

try:
   includes = env['ufc_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass

try:
   includes = env['amd_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass

try:
   lib_path = env['amd_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass

if useMPI:
    amd_libs = []
else:
   try:
      amd_libs = env['amd_libs']
   except KeyError:
      amd_libs = []

try:
   includes = env['blas_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass

try:
   lib_path = env['blas_lib_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass

try:
   blas_libs = env['blas_libs']
except KeyError:
   blas_libs = []

try:
   useNetCDF = env['useNetCDF']
except KeyError:
   useNetCDF = 'yes'
   pass

if not useNetCDF == 'yes':
   print "Warning: Installation is not configured with netCDF. Some I/O function may not be available."
   
if useNetCDF == 'yes': 
   try:
      includes = env['netCDF_path']
      env.Append(CPPPATH = [includes,])
   except KeyError:
      pass

   try:
      lib_path = env['netCDF_lib_path']
      if IS_WINDOWS_PLATFORM: env['ENV']['PATH']+=";"+lib_path
      env.Append(LIBPATH = [ lib_path, ])
   except KeyError:
      pass

   try:
      netCDF_libs = env['netCDF_libs']
   except KeyError:
      netCDF_libs = [ ]
else:
   netCDF_libs=[ ]

try:
   includes = env['boost_path']
   env.Append(CPPPATH = [includes,])
except KeyError:
   pass
try:
   lib_path = env['boost_libs_path']
   env.Append(LIBPATH = [lib_path,])
except KeyError:
   pass
try:
   boost_libs = env['boost_libs']
except KeyError:
   boost_libs = None
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
   python_libs = env['python_libs']
except KeyError:
   python_libs = None
try:
   doxygen_path = env['doxygen_path']
except KeyError:
   doxygen_path = None
try:
   epydoc_path = env['epydoc_path']
except KeyError:
   epydoc_path = None
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


try:
   src_zipfile = env.File(env['src_zipfile'])
except KeyError:
   src_zipfile = None
try:
   test_zipfile = env.File(env['test_zipfile'])
except KeyError:
   test_zipfile = None
try:
   examples_zipfile = env.File(env['examples_zipfile'])
except KeyError:
   examples_zipfile = None

try:
   src_tarfile = env.File(env['src_tarfile'])
except KeyError:
   src_tarfile = None
try:
   test_tarfile = env.File(env['test_tarfile'])
except KeyError:
   test_tarfile = None
try:
   examples_tarfile = env.File(env['examples_tarfile'])
except KeyError:
   examples_tarfile = None

try:
   guide_pdf = env.File(env['guide_pdf'])
except KeyError:
   guide_pdf = None

try:
   guide_html_index = env.File('index.htm',env['guide_html'])
except KeyError:
   guide_html_index = None

try:
   api_epydoc = env.Dir(env['api_epydoc'])
except KeyError:
   api_epydoc = None

try:
   api_doxygen = env.Dir(env['api_doxygen'])
except KeyError:
   api_doxygen = None

# Zipgets
env.Default(libinstall)
env.Default(incinstall)
env.Default(pyinstall)
### env.Default(exinstall) # ksteube this causes dependency error
env.Alias('release_src',[ src_zipfile, src_tarfile ])
env.Alias('release_tests',[ test_zipfile, test_tarfile])
env.Alias('release_examples',[ examples_zipfile, examples_tarfile])
env.Alias('examples_zipfile',examples_zipfile)
env.Alias('examples_tarfile',examples_tarfile)
env.Alias('api_epydoc',api_epydoc)
env.Alias('api_doxygen',api_doxygen)
env.Alias('guide_html_index',guide_html_index)
env.Alias('guide_pdf', guide_pdf)
env.Alias('docs',[ 'release_examples', 'guide_pdf', api_epydoc, api_doxygen, guide_html_index])
env.Alias('release', ['release_src', 'release_tests', 'docs'])
env.Alias('build_tests')    # target to build all C++ tests
env.Alias('build_py_tests') # target to build all python tests
env.Alias('build_all_tests', [ 'build_tests', 'build_py_tests' ] ) # target to build all python tests
env.Alias('run_tests', 'build_tests')   # target to run all C++ test
env.Alias('py_tests', 'build_py_tests') # taget to run all released python tests
env.Alias('all_tests', ['run_tests', 'py_tests']) # target to run all C++ and released python tests

# Python install - esys __init__.py
init_target = env.Command(pyinstall+'/__init__.py', None, Touch('$TARGET'))
env.Alias(init_target)

# Allow sconscripts to see the env
Export(["IS_WINDOWS_PLATFORM", "env", "incinstall", "libinstall", "pyinstall", "exinstall", "dodebug", "mkl_libs", "scsl_libs", "umf_libs", "amd_libs", "blas_libs", "netCDF_libs", "useNetCDF",
	"boost_libs", "python_libs", "doxygen_path", "epydoc_path", "papi_libs",
        "sys_libs", "test_zipfile", "src_zipfile", "test_tarfile", "src_tarfile", "examples_tarfile", "examples_zipfile",
        "guide_pdf", "guide_html_index", "api_epydoc", "api_doxygen", "useMPI" ])

# End initialisation section
# Begin configuration section
# adds this file and the scons option directore to the source tar
release_srcfiles=[env.File('SConstruct'),]+[ env.File(x) for x in glob.glob('scons/*.py') ]
release_testfiles=[env.File('README_TESTS'),]
env.Zip(src_zipfile, release_srcfiles)
env.Zip(test_zipfile, release_testfiles)
try:
   env.Tar(src_tarfile, release_srcfiles)
   env.Tar(test_tarfile, release_testfiles)
except AttributeError:
   pass
# Insert new components to be build here
# FIXME: might be nice to replace this verbosity with a list of targets and some
# FIXME: nifty python to create the lengthy but very similar env.Sconscript lines
# Third Party libraries
env.SConscript(dirs = ['tools/CppUnitTest/src'], build_dir='build/$PLATFORM/tools/CppUnitTest', duplicate=0)
# C/C++ Libraries
env.SConscript(dirs = ['paso/src'], build_dir='build/$PLATFORM/paso', duplicate=0)
# bruce is removed for now as it doesn't really do anything
# env.SConscript(dirs = ['bruce/src'], build_dir='build/$PLATFORM/bruce', duplicate=0)
env.SConscript(dirs = ['escript/src'], build_dir='build/$PLATFORM/escript', duplicate=0)
env.SConscript(dirs = ['esysUtils/src'], build_dir='build/$PLATFORM/esysUtils', duplicate=0)
env.SConscript(dirs = ['finley/src'], build_dir='build/$PLATFORM/finley', duplicate=0)
env.SConscript(dirs = ['modellib/py_src'], build_dir='build/$PLATFORM/modellib', duplicate=0)
env.SConscript(dirs = ['doc'], build_dir='build/$PLATFORM/doc', duplicate=0)
env.SConscript(dirs = ['pyvisi/py_src'], build_dir='build/$PLATFORM/pyvisi', duplicate=0)
env.SConscript(dirs = ['pycad/py_src'], build_dir='build/$PLATFORM/pycad', duplicate=0)

# added by Ben Cumming
env.SConscript(dirs = ['pythonMPI/src'], build_dir='build/$PLATFORM/pythonMPI', duplicate=0)
#env.SConscript(dirs = ['../test'], build_dir='../test/build', duplicate=0)
