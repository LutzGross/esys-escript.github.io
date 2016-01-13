# this code should be called by all esys13 scons builder and testing scripts
# it sets default values for relevant variables and overloads them with the 
# standard gcc settings by importing gcc_options. If there is 
# a <hostname>_options file it is imported to set platform specific settings

# this is the general set up for the esys scons system:
libinstall = None
pyinstall = None
incinstall = None
dodebug = 0

# locations of libs etc used by mkl
mkl_path = ''
mkl_lib_path = ''
mkl_libs = []

# locations of libs etc used by SCSL
scsl_path = ''
scsl_lib_path = ''
scsl_libs = []

# locations of libs etc used by UMFPACK
umfpack_path = ''
umfpack_lib_path = ''
umfpack_libs = []

# locations of include files for python
python_path = Dir('/usr/include')
python_lib_path = Dir('/usr/lib')
python_lib = [ 'python2.3' ]

# locations of PAPI
papi_path = ''
papi_lib_path = ''
papi_libs = []

# locations of libraries for boost
boost_path = '/usr/include'
boost_lib_path = '/usr/lib'
boost_lib = 'boost_python'

# locations of doc building executables
doxygen_path = ''
epydoc_path = ''
epydoc_pythonpath = ''

# names of c and c++ compilers to use
cc = 'gcc'
cxx = 'g++'

# c flags to use
cc_flags  = '-O0 -std=c99 -fpic -W -Wall -Wno-unknown-pragmas'
cc_flags_debug  = '-g -O0 -std=c99 -fpic -W -Wall -Wno-unknown-pragmas'

# c++ flags to use
cxx_flags  = '-O0 -ansi -fpic -W -Wall -Wno-unknown-pragmas'
cxx_flags_debug  = '-g -O0 -ansi -fpic -W -Wall -Wno-unknown-pragmas -DDOASSERT -DDOPROF'

# static library archiver flags to use
ar_flags = 'crus'

# system specific libraries to link with
sys_libs = []

#==== end of setting options ===========================================

import sys

# set esysroot 
options_dir = Dir(esysroot + '/scons')
if sys.path.count(str(options_dir))==0: sys.path.append(str(options_dir))

#
# ensure correct versions of python and scons
EnsurePythonVersion(2,3)
EnsureSConsVersion(0,96)

#
# retreive command-line arguments if any and overwrite settings in <hostname>_options
usegcc = 0
options = None
if ARGUMENTS.get('libinstall',0): libinstall = ARGUMENTS.get('libinstall',0)
if ARGUMENTS.get('pyinstall',0): pyinstall = ARGUMENTS.get('pyinstall',0)
if ARGUMENTS.get('incinstall',0): pyinstall = ARGUMENTS.get('incinstall',0)
if ARGUMENTS.get('debug',0): dodebug = 1
if ARGUMENTS.get('options',0): options = ARGUMENTS.get('options',0)
if ARGUMENTS.get('usegcc',0): usegcc = 1
#
# try to import <hostname>_options
try: 
    exec "from gcc_options import *"
except ImportError:
    pass

#
# try to import <hostname>_options
if usegcc==0:
   import socket
   from string import ascii_letters,digits
   hostname=""
   for s in socket.gethostname().split('.')[0]:
        if s in ascii_letters+digits:
           hostname+=s
        else:
           hostname+="_"
   try: 
      exec "from "+hostname+"_options import *"
   except ImportError:
      pass

#
# import additional options:
if options!=None:
  exec "from " + options + " import *"

#
# use debug options:
if dodebug==1:
     cxx_flags=cxx_flags_debug
     cc_flags=cc_flags_debug

#
# export configuration variables
#
Export(["esysroot"])
Export(["libinstall"])
Export(["incinstall"])
Export(["pyinstall"])
