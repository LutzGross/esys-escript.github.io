
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


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os, re, sys
try:
    import sysconfig
    USE_DISTUTILS = False
except ImportError:
    from distutils import sysconfig
    USE_DISTUTILS = True
from subprocess import PIPE, Popen
from SCons.Script.SConscript import Configure
from site_init import findLibWithHeader, detectModule
import subprocess

REQUIRED_BOOST = (1, 46)

def CheckComplexAcos(context):
    context.Message('Checking for working complex std::acos()... ')
    result = context.TryRun("""
#include <complex>
int main() { std::complex<double> x(0,3.14159265359), y(1.5707963,-1.8622957);
return std::abs(std::acos(x)-y) < 1e-6 ? 0:-1;}
""", '.cpp')
    # scons < 2.4 fix:
    if type(result)==tuple:
        result = result[0]
    context.Result(result)
    return result

def checkCompiler(env):
    conf = Configure(env.Clone(), custom_tests = {'CheckComplexAcos': CheckComplexAcos})
    if 'CheckCXX' in dir(conf): # exists since scons 1.1.0
        if not conf.CheckCXX():
            print("Cannot run C++ compiler '%s' (check config.log)" % (env['CXX']))
            env.Exit(1)
    else:
        if not conf.CheckFunc('printf', language='c++'):
            print("Cannot run C++ compiler '%s' (check config.log)" % (env['CXX']))
            env.Exit(1)

    conf.env['buildvars']['cxx']=conf.env['CXX']

    if conf.CheckFunc('gethostname', language='c++'):
        conf.env.Append(CPPDEFINES = ['HAVE_GETHOSTNAME'])

    if conf.CheckCXXHeader('byteswap.h'):
        checkhdr="""#include <byteswap.h>
#define SCbswap32() {int x=0;bswap_32(x);}"""
        if conf.CheckFunc('SCbswap32', header=checkhdr, language='c++'):
            conf.env.Append(CPPDEFINES = ['HAVE_BYTESWAP_H'])
    if conf.CheckCXXHeader('sys/endian.h'):
        conf.env.Append(CPPDEFINES = ['HAVE_SYS_ENDIAN_H'])
    if conf.CheckCXXHeader('libkern/OSByteOrder.h'):
        conf.env.Append(CPPDEFINES = ['HAVE_OSBYTEORDER_H'])

    if not conf.CheckComplexAcos():
        conf.env.Append(CPPDEFINES = ['ESYS_USE_BOOST_ACOS'])

    return conf.Finish()

def get_external_python_sympy(env,bin):
    if env['sympy']:
        import subprocess
        cmd=''
        cmd+='import sympy\n'
        cmd+='print(sympy.__version__)\n'
        sp=subprocess.Popen([bin, '-c', cmd], stdin=None, stderr=None, stdout=subprocess.PIPE)
        spVer=sp.stdout.readline()
        if hasattr(spVer, 'decode'):
            spVer=spVer.decode()
        import sys
        if sys.version_info[0] >= 3:
            spVer = spVer.strip()
        else:
            spVer = spVer.strip().split('.')
        quit=False
        ver1=''
        ver2=''
        count=0;
        for i in range(0,len(spVer)):
            x=spVer[i]
            if x.isdigit() is True:
                if count == 0:
                    ver1=ver1+spVer[i]
                else:
                    ver2=ver2+spVer[i]
            else:
                count=count+1
                if quit is True:
                    break
                else:
                    quit=True
                    continue
        version1=float(ver1)
        version2=float(ver2)
        if version1 == 0 and version2 < 7:
            env['sympy']=False
            env['warnings'].append("sympy version too old. Symbolic toolbox and nonlinear PDEs will not be available.")
            env.Append(CPPDEFINES = ['ESYS_NO_SYMPY'])
        #elif version1 == 1 and version2 > 2:
                                              #    env['sympy']=False
        #    env['warnings'].append("escript does not support sympy version 1.2 and higher. Found %s" % spVer)
        #    env.Append(CPPDEFINES = ['ESYS_NO_SYMPY'])
        else:
            env['sympy']=True
            env['warnings'].append("Found sympy version %s" % spVer)
        return env
    else:
        return env

def call_python_config(bin=None):
    import subprocess
    cmd=''
    cmd+='import subprocess\n'
    cmd+='import os\n'
    cmd+='import sys\n'
    if USE_DISTUTILS:
        cmd+='from distutils import sysconfig\n'
    else:
        cmd += 'import sysconfig\n'
    cmd+='pyversion=sysconfig.get_python_version()\n'
    cmd+='try:\n'
    cmd+='  sp=subprocess.Popen(["python"+pyversion+"-config","--ldflags"], stdout=subprocess.PIPE)\n'
    cmd+='except:\n'
    cmd+='  pythonroot=sys.exec_prefix+"/bin/"\n'
    cmd+='  sp=subprocess.Popen([pythonroot+"python"+pyversion+"-config","--ldflags"], stdout=subprocess.PIPE)\n'
    cmd+='d=sp.stdout.readline()\n'
    cmd+='if hasattr(d, "decode"):\n'
    cmd+='    d=d.decode()\n'
    cmd+='d=d.split()\n'
    cmd+="libdirs=[z[2:] for z in d if z.startswith('-L')]\n"
    cmd+="libs=[z[2:] for z in d if z.startswith('-lpython')]\n"
    cmd+="if len(libs) == 0:\n"
    cmd+='   libs = [ "python"+pyversion]\n'
    cmd+="target=''\n"
    cmd+="libname=''\n"
    cmd+="for d in libdirs:\n"
    cmd+="  for f in libs:\n"
    cmd+="    s=os.path.join(d,'lib'+f+'.so')\n"
    cmd+="    try:\n"
    cmd+="      dummy=os.stat(s)\n"
    cmd+="      if target=='':\n"
    cmd+="        target=d\n"
    cmd+="        libname=f\n"
    cmd+="    except Exception:\n"
    cmd+="      pass\n"
    cmd+="    s=os.path.join(d,'lib'+f+'.dylib')\n"
    cmd+="    try:\n"
    cmd+="      dummy=os.stat(s)\n"
    cmd+="      if target=='':\n"
    cmd+="        target=d\n"
    cmd+="        libname=f\n"
    cmd+="    except Exception:\n"
    cmd+="      pass\n"
    if bin is None:
       exec(cmd)
       ver=str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
       if USE_DISTUTILS:
           return (target,libname,ver, sysconfig.get_python_inc())
       else:
           return (target, libname, ver, sysconfig.get_path('include'))
    # run an external python to get its library location
    # yes we are starting another python just to run python?-config
    cmd+="if hasattr(target,'decode'):\n"
    cmd+="   target=target.decode()\n"
    cmd+="   libname=libname.decode()\n"
    cmd+="print(target)\n"
    cmd+="print(libname)\n"
    cmd+="import sys\n"
    cmd+="print(str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2]))\n"
    if USE_DISTUTILS:
        cmd += "print(sysconfig.get_python_inc())\n"
    else:
        cmd += "print(sysconfig.get_path('include'))\n"
    sp=subprocess.Popen([bin, '-c', cmd], stdin=None, stderr=None, stdout=subprocess.PIPE)
    target = sp.stdout.readline()
    libname = sp.stdout.readline()
    ver = sp.stdout.readline()
    pinc = sp.stdout.readline()
    if hasattr(target, "decode"):
        target, libname, ver, pinc=target.decode().strip(), libname.decode().strip(), ver.decode().strip(), pinc.decode().strip()
    else:
        target, libname, ver, pinc=target.strip(), libname.strip(), ver.strip(), pinc.strip()
    return (target, libname, ver, pinc)

def checkPython(env):
    # First we check to see if the config file has specified
    # where to find the file. Ideally, this should be automatic
    # but we need to deal with the case where python is not in its INSTALL
    # directory.
    # Use the python scons is running
    if env['pythoncmd'] == sys.executable:
        if env['IS_WINDOWS']:
            if USE_DISTUTILS:
                python_inc_path = sysconfig.get_python_inc()
            else:
                python_inc_path = sysconfig.get_path('include')
            python_lib_path=os.path.join(sysconfig.get_config_var('prefix'), 'libs')
            python_libs=['python%s%s'%(sys.version_info[0], sys.version_info[1])]
            verstring=".".join([str(i) for i in sys.version_info[:3]])
        else:
            (python_lib_path, python_libs,verstring, python_inc_path)=call_python_config(env['pythoncmd'])

    # if we want to use a python other than the one scons is running
    # Note: we assume scons is running python 2 in the following.
    else:
        if env['IS_WINDOWS']:
            cmd = "import os, sys\n"
            if USE_DISTUTILS:
                cmd += 'from distutils import sysconfig\n'
            else:
                cmd += 'import sysconfig\n'
            cmd += "print(sysconfig.get_python_inc())\n"
            cmd += "print(os.path.join(sysconfig.get_config_var('prefix'), 'libs'))\n"
            cmd += "print('python%s%s'%(sys.version_info[0], sys.version_info[1]))\n"
            cmd += "print('.'.join([str(i) for i in sys.version_info[:3]]))"
            pout = subprocess.Popen([env['pythoncmd'], '-c', cmd], stdout=subprocess.PIPE).stdout.read()
            if isinstance(pout, bytes):
                pout = pout.decode(sys.stdout.encoding)
            lines = pout.split('\n')
            python_inc_path = lines[0].strip()
            python_lib_path = lines[1].strip()
            python_libs = [lines[2].strip()]
            verstring = lines[3].strip()
        else:
            (python_lib_path, python_libs,verstring, python_inc_path)=call_python_config(env['pythoncmd'])
    if isinstance(python_inc_path, bytes):
        python_inc_path=python_inc_path.decode()
    if isinstance(python_lib_path, bytes):
        python_lib_path=python_lib_path.decode()
    env['python_version'] = verstring
    try:
        ispython3 = (verstring[0] == '3')
    except:
        ispython3 = True
    if ispython3:
        env.Append(CPPDEFINES=['ESPYTHON3'])
    env['buildvars']['python_version'] = verstring
    env['buildvars']['python'] = env['pythoncmd']
    # Check for an override from the config file.
    # Ideally, this should be automatic but we need to deal with the case
    # where python is not in its INSTALL directory
    if env['pythonlibpath'] :
        python_lib_path = env['pythonlibpath']

    if env['pythonincpath'] :
        python_inc_path = env['pythonincpath']

    if env['pythonlibname'] :
        python_libs = env['pythonlibname']

    conf = Configure(env.Clone())

    if env['sysheaderopt'] :
        conf.env.AppendUnique(CCFLAGS= [ (env['sysheaderopt'] , python_inc_path) ] )
    else:
        conf.env.AppendUnique(CPPPATH = [python_inc_path])

    conf.env.AppendUnique(LIBPATH = [python_lib_path])
    conf.env.AppendUnique(LIBS = python_libs)
    # The wrapper script needs to find the libs
    conf.env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], python_lib_path)

    if not conf.CheckCXXHeader('Python.h'):
        print("Cannot find python include files (tried 'Python.h' in directory %s)" % (python_inc_path))
        env.Exit(1)
    if not env['IS_WINDOWS']:
        if not conf.CheckFunc('Py_Exit', language='c++'):
            print("Cannot find python library method Py_Exit (tried %s in directory %s)" % (python_libs, python_lib_path))
            env.Exit(1)

    return conf.Finish()

def checkCudaVersion(env):
    # NVCC availability is already checked in the Tool file
    p = Popen([env['NVCC'], '-V'], stdout=PIPE)
    out,_ = p.communicate()
    env['nvcc_version'] = '(unknown version)'
    for line in out.split('\n'):
        if 'release' in line:
            version = line[line.find('release'):].strip()
            env['nvcc_version'] = version
            break
    env['buildvars']['nvcc']=env['NVCC']
    return env

def checkBoost(env):

    boost_inc_path, boost_lib_path = findLibWithHeader(env, env['boost_libs'], 'boost/python.hpp', env['boost_prefix'], lang='c++')

    if env['sysheaderopt'] == '':
        env.AppendUnique(CPPPATH = [ boost_inc_path ])
    else:
        # This is required because we can't -isystem /usr/include since it
        # breaks std includes
        if os.path.normpath(boost_inc_path) == '/usr/include':
            env.Append(CCFLAGS= [ (env['sysheaderopt'], os.path.join(boost_inc_path,'boost')) ] )
        else:
            env.Append(CCFLAGS= [ (env['sysheaderopt'], boost_inc_path ) ] )

    env.AppendUnique(LIBPATH = [ boost_lib_path ])
    env.AppendUnique(LIBS = env['boost_libs'])
    env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], boost_lib_path)

    # Try to extract the boost version from version.hpp
    boosthpp=open(os.path.join(boost_inc_path, 'boost', 'version.hpp'))
    boostversion='unknown'
    for line in boosthpp:
        ver=re.match(r'#define BOOST_VERSION (\d+)',line)
        if ver:
            boostversion=ver.group(1)
            boostversion = int(boostversion)
            maj = boostversion/100000
            minor = (boostversion/100)%1000
            sub = boostversion % 100
            env['boost_version'] = "%d.%d.%d"%(maj,minor,sub)
            if maj <= REQUIRED_BOOST[0] and minor < REQUIRED_BOOST[1]:
                print("The boost version referenced must be at least version %d.%d "%REQUIRED_BOOST + "(have %d.%d.%d)"%(maj,minor,sub))
                env.Exit(1)
    boosthpp.close()
    env['buildvars']['boost_inc_path']=boost_inc_path
    env['buildvars']['boost_lib_path']=boost_lib_path
    env['buildvars']['boostversion']=boostversion

    #=======================
    # Check for the boost numpy library
    env['have_boost_numpy']=False
    if boostversion >= 106300 and not env['disable_boost_numpy'] :
        pyv = env['python_version'].split(".")
        if env["IS_OSX"]:
            libname = f'boost_numpy{pyv[0]}{pyv[1]}'
        else:
            libname = f'boost_numpy{pyv[0]}{pyv[1]}'
        try:
            boost_numpy_inc_path, boost_numpy_lib_path = \
                findLibWithHeader(env, [ libname ] + env['boost_libs'], os.path.join('boost', 'python', 'numpy.hpp'), env['boost_prefix'], lang='c++')
            env.AppendUnique(LIBS=libname)
            env.AppendUnique(boost_libs=libname)
            env.AppendUnique(CPPPATH=boost_numpy_inc_path)
            env.AppendUnique(LIBPATH=boost_numpy_lib_path)
            env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], boost_numpy_lib_path)
            env.Append(CPPDEFINES=['ESYS_HAVE_BOOST_NUMPY'])
            env['have_boost_numpy'] = True
        except:
            print("Warning: Could not find boost/python/numpy.hpp. Building without numpy support.")

    if boostversion >= 107000:
        env.Append(CPPDEFINES=['ESYS_DEPRECATED_BOOST_ENDIAN'])
    if boostversion >= 107200:
        env.Append(CPPDEFINES=['ESYS_MOVED_BOOST_ENDIAN'])

    return env

def checkNumpy(env):
    if not detectModule(env, 'numpy'):
        print("Cannot import numpy. If it is installed try setting your PYTHONPATH and probably %s"%env['LD_LIBRARY_PATH_KEY'])
        env.Exit(1)

    ## check for numpy header (optional)
    conf = Configure(env.Clone())
    numpy_h = False
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
            # make a copy of CPPPATH so it can be restored if header check fails
            cpp_path_old = conf.env.get('CPPPATH', []).copy()
            conf.env.Append(CPPPATH = [conda_prefix+'/Lib/site-packages/numpy/core/include'])
            if conf.CheckCXXHeader(['Python.h','numpy/ndarrayobject.h']):
                numpy_h = True
            else:
                conf.env['CPPPATH'] = cpp_path_old
    elif env['IS_OSX']:
        # make a copy of CPPPATH so it can be restored if header check fails
        cpp_path_old = conf.env.get('CPPPATH', []).copy()
        pypth = env['pythoncmd'][:env['pythoncmd'].index('bin')-1]
        conf.env.Append(CPPPATH=[
            f'{pypth}/lib/python{env["python_version"][:env["python_version"].find(".",2)]}/site-packages/numpy/core/include'])
        if conf.CheckCXXHeader(['Python.h', 'numpy/ndarrayobject.h']):
            numpy_h = True
        else:
            conf.env['CPPPATH'] = cpp_path_old
    else:
        cpp_path_old = conf.env.get('CPPPATH', []).copy()
        if conf.CheckCXXHeader(['Python.h','numpy/ndarrayobject.h']):
            numpy_h = True
        else:
            conf.env.Append(CPPPATH=['/usr/lib/x86_64-linux-gnu/python3-numpy/numpy/_core/include'])
            if conf.CheckCXXHeader(['Python.h', 'numpy/ndarrayobject.h']):
                numpy_h = True
            else:
                # Extract Python prefix from pythoncmd
                if 'bin' in env['pythoncmd']:
                    pypth = env['pythoncmd'][:env['pythoncmd'].index('bin')-1]
                else:
                    # pythoncmd is just the command name, use sys.prefix
                    import sys
                    pypth = sys.prefix
                conf.env.Append(CPPPATH=[
                    f'{pypth}/lib/python{env["python_version"][:env["python_version"].find(".",2)]}/site-packages/numpy/_core/include'])
                if conf.CheckCXXHeader(['Python.h', 'numpy/ndarrayobject.h']):
                    numpy_h = True
                else:
                    conf.env['CPPPATH'] = cpp_path_old
    conf.env['numpy_h'] = numpy_h
    if numpy_h:
        conf.env.Append(CPPDEFINES = ['ESYS_HAVE_NUMPY_H'])

    return conf.Finish()

def checkCUDA(env):
    try:
        cuda_inc_path,cuda_lib_path=findLibWithHeader(env, 'cudart', 'thrust/version.h', env['cuda_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [cuda_inc_path])
        env.AppendUnique(LIBPATH = [cuda_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], cuda_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_CUDA'])
        env['cuda']=True
    except:
        env['cuda']=False
    return env

def checkCppUnit(env):
    try:
        cppunit_inc_path,cppunit_lib_path=findLibWithHeader(env, env['cppunit_libs'], 'cppunit/TestFixture.h', env['cppunit_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [cppunit_inc_path])
        env.AppendUnique(LIBPATH = [cppunit_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], cppunit_lib_path)
        env['cppunit']=True
    except:
        env['cppunit']=False
    return env

def checkOptionalModules(env):
    ######## scipy
    if not detectModule(env, 'scipy'):
        env['warnings'].append("Cannot import scipy.")

    ######## sympy
    if env['sympy']:
        if detectModule(env, 'sympy'):
            env['sympy'] = True
            env['warnings'].append("Found sympy.")
        else:
            env.Append(CPPDEFINES = ['ESYS_NO_SYMPY'])
            env['sympy']=False
            env['warnings'].append("Could not find sympy")
    else:
        env['sympy']=False
        env.Append(CPPDEFINES = ['ESYS_NO_SYMPY'])
    return env

def checkForTrilinos(env):
    if not env['trilinos']:
        env['trilinos_version']='none'
        return env
    havelibs = (len(env['trilinos_libs']) > 0)
    trilinos_inc_path, trilinos_lib_path=findLibWithHeader(env,
            env['trilinos_libs'], 'Amesos2.hpp',
            env['trilinos_prefix'], lang='c++', try_link=havelibs)

    # Try to extract the trilinos version from Trilinos_version.h
    versionh = open(os.path.join(trilinos_inc_path, 'Trilinos_version.h'))
    env['trilinos_version'] = 'unknown'
    for line in versionh:
        ver = re.match(r'#define TRILINOS_MAJOR_MINOR_VERSION (\d+)', line)
        if ver:
            trilinos_version = ver.group(1)
            trilinos_version = int(trilinos_version)
            major = int(str(trilinos_version)[:2])
            minor = int(str(trilinos_version)[2:4])
            subminor = int(str(trilinos_version)[4:6])
            env['trilinos_version'] = str(major) + "." + str(minor) + "." + str(subminor)
            if major < 14:
                raise RuntimeError('Trilinos version greater than 14 expected.')
            #            if major == 14 and minor <2 :
            #                env.Append(CPPDEFINES = ['ESYS_TRILINOS_14'])
            #            elif major == 14 and minor >=2:
            #                env.Append(CPPDEFINES = ['ESYS_TRILINOS_14_2'])
            #            else:
            env.Append(CPPDEFINES=[f'ESYS_TRILINOS_{major}'])

    packages = ['Tpetra', 'Kokkos', 'Belos', 'Amesos2', 'Ifpack2', 'MueLu']

    dependencies = ['Amesos2.hpp', 'Amesos2_Solver_decl.hpp', 'BelosSolverFactory.hpp', 'BelosSolverManager.hpp', \
                    'BelosTFQMRIter.hpp', 'BelosTFQMRSolMgr.hpp', 'BelosTpetraAdapter.hpp', 'BelosTypes.hpp', \
                    'Ifpack2_Factory.hpp',  'Tpetra_KokkosCompat_DefaultNode.hpp', \
                    'MatrixMarket_Tpetra.hpp', 'MueLu_CreateTpetraPreconditioner.hpp', \
                    'Teuchos_DefaultComm.hpp', 'Teuchos_ParameterList.hpp', 'Tpetra_BlockCrsMatrix.hpp', \
                    'Teuchos_Comm.hpp', 'Teuchos_TimeMonitor.hpp', 'Tpetra_CrsMatrix_decl.hpp', \
                    'Tpetra_BlockCrsMatrix_decl.hpp',  \
                    'Tpetra_CrsGraph.hpp', 'Tpetra_CrsMatrix.hpp', 'Tpetra_RowMatrix.hpp', \
                    'TpetraExt_TripleMatrixMultiply_def.hpp', 'Tpetra_BlockVector.hpp', \
                    'Tpetra_Vector.hpp', 'Trilinos_version.h', 'Tpetra_BlockCrsMatrix_Helpers.hpp']

    env.AppendUnique(LIBPATH = [trilinos_lib_path])
    env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], trilinos_lib_path)
    env['buildvars']['trilinos_lib_path']=trilinos_lib_path


    if major > 16 or ( major == 16 and minor > 0) :
        dependencies +=[ os.path.join('kokkos', 'Kokkos_Core.hpp' ) ]
        env.AppendUnique(CPPPATH=[trilinos_inc_path, os.path.join(trilinos_inc_path, 'kokkos')])
        env['buildvars']['trilinos_inc_path'] = [trilinos_inc_path, os.path.join(trilinos_inc_path, 'kokkos')]
    else:
        env.AppendUnique(CPPPATH=[trilinos_inc_path])
        env['buildvars']['trilinos_inc_path'] = [ trilinos_inc_path ]
        dependencies +=['Kokkos_Core.hpp']

    print("Looking for the Trilinos headers...")
    for check in dependencies:
        print("Checking for %s... %s" % (check, "yes" if os.path.isfile(os.path.join(trilinos_inc_path,check)) else "no") )
        if not os.path.isfile(os.path.join(trilinos_inc_path,check)):
            print("Could not find the  Trilinos header file %s (tried looking in directory %s)" % (check, trilinos_inc_path))
            env.Exit(1)
    # Check for required Tpetra headers (standard versions required for Trilinos 16.2+)
    required_tpetra_headers = [
        'Tpetra_BlockCrsMatrix.hpp',
        'Tpetra_BlockCrsMatrix_Helpers.hpp',
        'Tpetra_BlockVector.hpp'
    ]
    for header in required_tpetra_headers:
        header_path = os.path.join(trilinos_inc_path, header)
        if os.path.isfile(header_path):
            print("Checking for %s... yes" % header)
        else:
            raise RuntimeError('Could not locate the Trilinos header %s (required for Trilinos 16.2+)' % header)

    library_list=['amesos2', 'amesos', 'anasaziepetra', 'anasazi', 'anasazitpetra', 'aztecoo', \
                  'belosepetra', 'belos', 'belostpetra', 'belosxpetra', 'epetraext', 'epetra', \
                  'galeri', 'galeri', 'ifpack2', 'ifpack2', 'ifpack', 'intrepid2', 'isorropia',\
                  'kokkosalgorithms', 'kokkoscontainers', 'kokkoscore', 'kokkoskernels', 'kokkossimd', \
                  'kokkostsqr', 'komplex', 'ml', 'ModeLaplace', 'muelu', 'muelu', 'pamgen_extras',\
                  'pamgen', 'rtop', 'sacado', 'shards', 'shylu_nodehts', 'simpi', 'stratimikosamesos2',\
                  'stratimikosamesos', 'stratimikosaztecoo', 'stratimikosbelos', 'stratimikosifpack', \
                  'stratimikosml', 'stratimikos', 'tacho', 'teko', 'teuchoscomm', 'teuchoscore', \
                  'teuchoskokkoscomm', 'teuchoskokkoscompat', 'teuchosnumerics', 'teuchosparameterlist', \
                  'teuchosparser', 'teuchosremainder', 'thyracore', 'thyraepetraext', 'thyraepetra', \
                  'thyratpetra', 'tpetraclassic', 'tpetraext', 'tpetrainout', 'tpetra', 'xpetra', 'xpetra-sup',
                  'trilinosss', 'triutils', 'xpetra', 'xpetra', 'zoltan', 'zoltan2']

    if not havelibs:
        libs = []
        for file in os.listdir(trilinos_lib_path):
            if file.endswith(".so"):
                for x in range(0, len(library_list)):
                    if file[3:-3] == library_list[x] or file[12:-3] == library_list[x]:
                        libs.append(file[3:-3])
            elif file.endswith(".dylib"):
                for x in range(0, len(library_list)):
                    if file[3:-6] == library_list[x] or file[12:-6] == library_list[x]:
                        libs.append(file[3:-6])
        env['trilinos_libs'] = libs
    # Trilinos version 14 and higher
    # if major >= 14:
    #     if not havelibs:
    #         libs = []
    #         for file in os.listdir(trilinos_lib_path):
    #             if file[-3:] == ".so":
    #                 libs.append(file[3:-3])
    #             if file[-3:] == ".dylib": # macOS
    #                 libs.append(file[3:-6])
    #         env['trilinos_libs'] = libs
    # else:
    #     if not havelibs:
    #         libs = []
    #         for pk in packages:
    #             # find out what libraries to link with...
    #             makefile = os.path.join(trilinos_inc_path, 'Makefile.export.%s'%pk)
    #             try:
    #                 for l in open(makefile, 'r').readlines():
    #                     if l.startswith("%s_LIBRARIES"%pk): # or l.startswith("Trilinos_TPL_LIBRARIES"):
    #                         lst = l.split('=')[1].strip().split()
    #                         lst = [e.replace('-l','',1) for e in lst]
    #                         libs += lst
    #                     elif l.startswith("%s_TPL_INCLUDE_DIRS"%pk):
    #                         lst = l.split('=')[1].strip().split()
    #                         lst = [e.replace('-I','',1) for e in lst]
    #                         env.AppendUnique(CPPPATH = lst)

    #             except Exception as e:
    #                 raise RuntimeError('Error reading Trilinos export Makefile\n%s'%(e))
    #         env['trilinos_libs'] = libs

    pytri_path=trilinos_lib_path+'/python'+str(sys.version_info[0])+'.'+str(sys.version_info[1])+'/site-packages/PyTrilinos'
    pytri_test_file=os.path.join(pytri_path,'Tpetra.pyc')
    if os.path.isfile(pytri_test_file):
        env.Append(PYTHONPATH=pytri_path)
    paths=sys.path
    for i in range(0,paths.__len__()):
        env.Append(PYTHONPATH=paths[i])
    env.Append(CPPDEFINES = ['ESYS_HAVE_TRILINOS'])
    # env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], trilinos_lib_path)
    # env['buildvars']['trilinos_inc_path']=trilinos_inc_path
    # env['buildvars']['trilinos_lib_path']=trilinos_lib_path
    env['buildvars']['trilinos']=int(env['trilinos'])
    return env

def checkOptionalLibraries(env):
    env['usempi'] = env['mpi']!='none'
    ######## hdf5
    hdf5_inc_path = ''
    hdf5_lib_path = ''
    if env['hdf5']:
        env.Append(CPPDEFINES=['HDF5'])
        hdf5_inc_path, hdf5_lib_path = findLibWithHeader(env, env['hdf5_libs'], 'H5Cpp.h',
                                                                 env['hdf5_prefix'], lang='c++')

        env.AppendUnique(CPPPATH=[hdf5_inc_path])
        env.AppendUnique(LIBPATH=[hdf5_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], hdf5_lib_path)
        env.Append(CPPDEFINES=['ESYS_HAVE_HDF5'])
        env['buildvars']['hdf5_inc_path'] = hdf5_inc_path
        env['buildvars']['hdf5_lib_path'] = hdf5_lib_path
    env['buildvars']['hdf5'] = int(env['hdf5'])
    ######## MKL
    mkl_inc_path=''
    mkl_lib_path=''
    if env['mkl']:
        mkl_inc_path,mkl_lib_path=findLibWithHeader(env, env['mkl_libs'], 'mkl_pardiso.h', env['mkl_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [mkl_inc_path])
        env.AppendUnique(LIBPATH = [mkl_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], mkl_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_MKL'])
        env['buildvars']['mkl_inc_path']=mkl_inc_path
        env['buildvars']['mkl_lib_path']=mkl_lib_path
    env['buildvars']['mkl']=int(env['mkl'])

    ######## UMFPACK
    umfpack_inc_path=''
    umfpack_lib_path=''
    if env['umfpack']:
        umfpack_inc_path,umfpack_lib_path=findLibWithHeader(env, env['umfpack_libs'], 'umfpack.h', env['umfpack_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [umfpack_inc_path])
        env.AppendUnique(LIBPATH = [umfpack_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], umfpack_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_UMFPACK'])
        env['buildvars']['umfpack_inc_path']=umfpack_inc_path
        env['buildvars']['umfpack_lib_path']=umfpack_lib_path
    env['buildvars']['umfpack']=int(env['umfpack'])

    ######## Sequential MUMPS (works with MPI builds)
    mumps_seq_inc_path=''
    mumps_seq_lib_path=''
    if env['mumps_seq']:
        # Sequential MUMPS uses mumps_seq/mpi.h
        mumps_seq_inc_path, mumps_seq_lib_path = findLibWithHeader(env, env['mumps_seq_libs'], 'mumps_seq/mpi.h', env['mumps_seq_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [mumps_seq_inc_path])
        env.AppendUnique(LIBPATH = [mumps_seq_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], mumps_seq_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_MUMPS'])
        env['buildvars']['mumps_seq_inc_path']=mumps_seq_inc_path
        env['buildvars']['mumps_seq_lib_path']=mumps_seq_lib_path
    env['buildvars']['mumps_seq']=int(env['mumps_seq'])

    ######## LAPACK
    lapack_inc_path=''
    lapack_lib_path=''
    flavour = 'none'
    env['uselapack'] = False
    if env['lapack'] != 0:
        # not explicitly disabled so run the checks
        if env['longindices']:
            # you want longindices + lapack? sorry.
            if env['lapack'] == 1:
                print("LAPACK requires index type = int. Set longindices to False or disable LAPACK.")
                env.Exit(1)
        else:
            if env['mkl']:
                # we detected MKL so try the MKL header+libs
                flavour = 'mkl'
                header = 'mkl_lapack.h'
                prefix = env['mkl_prefix']
                if len(env['lapack_libs']) == 0:
                    libs = env['mkl_libs']
                else:
                    libs = env['lapack_libs']
            else:
                # try for lapacke
                flavour = 'lapacke'
                header = 'lapacke.h'
                prefix = env['lapack_prefix']
                if len(env['lapack_libs']) == 0:
                    libs = ['lapacke']
                else:
                    libs = env['lapack_libs']

            try:
                lapack_inc_path,lapack_lib_path=findLibWithHeader(env, libs, header, prefix, lang='c++')
                env['lapack_libs'] = libs
                env['uselapack'] = True
                env.AppendUnique(CPPPATH = [lapack_inc_path])
                env.AppendUnique(LIBPATH = [lapack_lib_path])
                env.Append(CPPDEFINES = ['ESYS_HAVE_LAPACK'])
                if flavour == 'mkl':
                    env.AppendUnique(CPPDEFINES = ['ESYS_MKL_LAPACK'])
                env['buildvars']['lapack_inc_path']=lapack_inc_path
                env['buildvars']['lapack_lib_path']=lapack_lib_path
            except:
                if env['lapack'] == 1:
                    raise
                # lapack was set to auto-detect so not a fatal error
                flavour = 'none'

    env['lapack'] = flavour
    env['buildvars']['lapack'] = flavour

    ######## Silo
    silo_inc_path=''
    silo_lib_path=''
    if env['silo']:
        silo_inc_path,silo_lib_path=findLibWithHeader(env, env['silo_libs'], 'silo.h', env['silo_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [silo_inc_path])
        env.AppendUnique(LIBPATH = [silo_lib_path])
        env.Append(CPPDEFINES = ['ESYS_HAVE_SILO'])
        env['buildvars']['silo_inc_path']=silo_inc_path
        env['buildvars']['silo_lib_path']=silo_lib_path
    env['buildvars']['silo']=int(env['silo'])

    ######## NetCDF
    netcdf_inc_path=''
    netcdf_lib_path=''
    if env['netcdf']:
        netcdf_inc_path,netcdf_lib_path=findLibWithHeader(env, env['netcdf_libs'], 'ncFile.h', env['netcdf_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [netcdf_inc_path])
        env.AppendUnique(LIBPATH = [netcdf_lib_path])
        env.Append(CPPDEFINES = ['ESYS_HAVE_NETCDF4'])
        env['buildvars']['netcdf_inc_path']=netcdf_inc_path
        env['buildvars']['netcdf_lib_path']=netcdf_lib_path
    env['buildvars']['netcdf']=int(bool(env['netcdf']))

    ######## VisIt
    visit_inc_path=''
    visit_lib_path=''
    if env['visit']:
        visit_inc_path,visit_lib_path=findLibWithHeader(env, env['visit_libs'], 'VisItControlInterface_V2.h', env['visit_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [visit_inc_path])
        env.AppendUnique(LIBPATH = [visit_lib_path])
        env['buildvars']['visit_inc_path']=visit_inc_path
        env['buildvars']['visit_lib_path']=visit_lib_path
    env['buildvars']['visit']=int(env['visit'])

    ######## MPI

    mpi_inc_path=''
    mpi_lib_path=''
    if env['usempi']:
        mpi_inc_path,mpi_lib_path=findLibWithHeader(env, env['mpi_libs'], 'mpi.h', env['mpi_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [mpi_inc_path])
        env.AppendUnique(LIBPATH = [mpi_lib_path])
        env.AppendUnique(LIBS = env['mpi_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], mpi_lib_path)
        env.Append(CPPDEFINES = ['ESYS_MPI', 'MPI_NO_CPPBIND', 'MPICH_IGNORE_CXX_SEEK'])

        if env['mpi'] == 'OPENMPI':
            # try to get version for correct launcher arguments
            try:
                p = Popen(['orterun', '-V'], stdout=PIPE, stderr=PIPE)
                o,e = p.communicate()
                try:
                    e=e.decode()
                    ver = e.split('\n')[0].split()[-1]
                except IndexError:
                    o=o.decode()
                    ver = o.split('\n')[0].split()[-1]
                if len(ver) > 0:
                    env['orte_version'] = ver
            except OSError:
                pass

        env['buildvars']['mpi_inc_path']=mpi_inc_path
        env['buildvars']['mpi_lib_path']=mpi_lib_path
    env['buildvars']['mpi']=env['mpi']

    ######## METIS
    metis_inc_path=''
    metis_lib_path=''
    if env['metis']:
        metis_inc_path,metis_lib_path=findLibWithHeader(env, env['metis_libs'], 'metis.h', env['metis_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [metis_inc_path])
        env.AppendUnique(LIBPATH = [metis_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], metis_lib_path)

        # Try to extract the metis version from metis.h
        header=open(os.path.join(metis_inc_path, 'metis.h')).readlines()
        major,minor,sub = None,None,None
        for line in header:
            ver=re.match(r'#define METIS_VER_MAJOR\s*(\d+)',line)
            if ver:
                major = int(ver.group(1))
                continue
            ver=re.match(r'#define METIS_VER_MINOR\s*(\d+)',line)
            if ver:
                minor = int(ver.group(1))
                continue
            ver=re.match(r'#define METIS_VER_SUBMINOR\s*(\d+)',line)
            if ver:
                sub = int(ver.group(1))
                continue
        if major is not None:
            env['metis_version'] = "%d.%d.%d"%(major,minor if minor else 0,sub if sub else 0)
        else:
            env['metis_version'] = "unknown"
    else:
        env['metis_version'] = "N/A"
    env['buildvars']['metis']=int(env['metis'])
    env['buildvars']['metis_inc_path']=metis_inc_path
    env['buildvars']['metis_lib_path']=metis_lib_path

    ######## ParMETIS
    if not env['usempi']: env['parmetis'] = False
    parmetis_inc_path=''
    parmetis_lib_path=''
    if env['parmetis']:
        parmetis_inc_path,parmetis_lib_path=findLibWithHeader(env, env['parmetis_libs'], 'parmetis.h', env['parmetis_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [parmetis_inc_path])
        env.AppendUnique(LIBPATH = [parmetis_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], parmetis_lib_path)

        # Try to extract the parmetis version from parmetis.h
        header=open(os.path.join(parmetis_inc_path, 'parmetis.h')).readlines()
        major,minor,sub = None,None,None
        for line in header:
            ver=re.match(r'#define PARMETIS_MAJOR_VERSION\s*(\d+)',line)
            if ver:
                major = int(ver.group(1))
                continue
            ver=re.match(r'#define PARMETIS_MINOR_VERSION\s*(\d+)',line)
            if ver:
                minor = int(ver.group(1))
                continue
            ver=re.match(r'#define PARMETIS_SUBMINOR_VERSION\s*(\d+)',line)
            if ver:
                sub = int(ver.group(1))
                continue
        if major is not None:
            env['parmetis_version'] = "%d.%d.%d"%(major,minor,0 if sub is None else sub)
            if env['longindices']:
                # ParMETIS version 3.x does not support 64-bit indices
                if major < 4:
                    print("Sorry, cannot use ParMETIS version < 4.0 with 64-bit index types. Set longindices to False or disable ParMETIS.")
                    env.Exit(1)
                else:
                    # check if ParMETIS was built with 64-bit indices
                    conf = Configure(env.Clone())
                    idxsize=conf.CheckTypeSize('idx_t', '#include <parmetis.h>', 'C++')
                    if idxsize != 8:
                        print("Sorry, ParMETIS was not compiled with 64-bit indices. Set longindices to False or disable/rebuild ParMETIS.")
                        env.Exit(1)
        else:
            env['parmetis_version'] = "unknown"

        env.Append(CPPDEFINES = ['ESYS_HAVE_PARMETIS'])
        env['buildvars']['parmetis_inc_path']=parmetis_inc_path
        env['buildvars']['parmetis_lib_path']=parmetis_lib_path
    env['buildvars']['parmetis']=int(env['parmetis'])


    ######## zlib
    if env['zlib']:
        try:
            zlib_inc_path, zlib_lib_path = findLibWithHeader(env, env['zlib_libs'], 'zlib.h', env['zlib_prefix'], lang='c')
            env.Append(CPPDEFINES = ['ESYS_HAVE_ZLIB'])
            env.AppendUnique(CPPPATH=[zlib_inc_path])
            env.AppendUnique(LIBPATH=[zlib_lib_path])
            env.AppendUnique(LIBS=env['zlib_libs'])
            env['buildvars']['zlib_inc_path'] = zlib_inc_path
            env['buildvars']['zlib_lib_path'] = zlib_lib_path
        except RuntimeError as e:
            env['zlib'] = False
    env['buildvars']['zlib']=int(env['zlib'])
    env['buildvars']['sympy']=int(env['sympy'])
    ######## boost::iostreams
    if env['compressed_files']:
        try:
            boost_inc_path, boost_lib_path = findLibWithHeader(env, env['compression_libs'], 'boost/iostreams/filter/gzip.hpp', env['boost_prefix'], lang='c++')
            env.Append(CPPDEFINES = ['ESYS_HAVE_BOOST_IO'])
        except RuntimeError as e:
            env['compressed_files'] = False
    env['buildvars']['compressed_files']=int(env['compressed_files'])

    ######## Trilinos
    env = checkForTrilinos(env)

    ######## mpi4py
    if env['mpi4py']:
        try:
            import mpi4py
            mpi4py_inc_path=mpi4py.get_include()
            env.AppendUnique(CPPPATH = [mpi4py_inc_path])
            env.Append(CPPDEFINES = ['ESYS_HAVE_MPI4PY'])
            env['mpi4py'] = True
        except RuntimeError as e:
            env['warnings'].append("Could not import mpi4py.")

    return env

def checkPDFLatex(env):
    if 'PDF' in dir(env) and '.tex' in env.PDF.builder.src_suffixes(env):
        env['pdflatex']=True
    else:
        env['pdflatex']=False
    return env
