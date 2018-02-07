
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import os, re, sys
from distutils import sysconfig
from subprocess import PIPE, Popen
from SCons.Script.SConscript import Configure
from site_init import findLibWithHeader, detectModule

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

def get_external_python_sympy(bin):
    import subprocess
    cmd=''
    cmd+='import sympy\n'
    cmd+='print(sympy.__version__)\n'
    sp=subprocess.Popen([bin, '-c', cmd], stdin=None, stderr=None, stdout=subprocess.PIPE)
    ver=sp.stdout.readline().strip().split('.')
    if int(ver[0]) == 0 and int(ver[1]) < 7:
        env['sympy'] = False
        env['warnings'].append("sympy version is too old.")

def call_python_config(bin=None):
    import subprocess
    cmd=''
    cmd+='import subprocess\n'
    cmd+='import os\n'
    cmd+='from distutils import sysconfig\n'
    cmd+='pyversion=sysconfig.get_python_version()\n'
    cmd+='sp=subprocess.Popen(["python"+pyversion+"-config","--ldflags"], stdout=subprocess.PIPE)\n'
    cmd+='d=sp.stdout.readline().split()\n'
    cmd+="libdirs=[z[2:] for z in d if z.startswith(b'-L')]\n"
    cmd+="libs=[z[2:] for z in d if z.startswith(b'-lpython')]\n"
    cmd+="target=''\n"
    cmd+="libname=''\n"
    cmd+="for d in libdirs:\n"
    cmd+="  for f in libs:\n"
    cmd+="    s=os.path.join(d,b'lib'+f+b'.so')\n"
    cmd+="    try:\n"
    cmd+="      dummy=os.stat(s)\n"
    cmd+="      if target=='':\n"
    cmd+="        target=d\n"
    cmd+="        libname=f\n"
    cmd+="    except Exception:\n"
    cmd+="      pass\n"
    cmd+="    s=os.path.join(d,b'lib'+f+b'.dylib')\n"
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
       return (target,libname,ver, sysconfig.get_python_inc())
    # run an external python to get its library location
    # yes we are starting another python just to run python?-config
    cmd+="if 'decode' in dir(target):\n"
    cmd+="   target=target.decode()\n"
    cmd+="   libname=libname.decode()\n"
    cmd+="print(target)\n"
    cmd+="print(libname)\n"
    cmd+="import sys\n"
    cmd+="print(str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2]))\n"
    cmd+="print(sysconfig.get_python_inc())\n"
    sp=subprocess.Popen([bin, '-c', cmd], stdin=None, stderr=None, stdout=subprocess.PIPE)
    target=sp.stdout.readline().strip()
    libname=sp.stdout.readline().strip()
    ver=sp.stdout.readline().strip()
    pinc=sp.stdout.readline().strip()
    return (target, libname, ver, pinc)

def checkPython(env):
    # First we check to see if the config file has specified
    # where to find the file. Ideally, this should be automatic
    # but we need to deal with the case where python is not in its INSTALL
    # directory.
    # Use the python scons is running
    if env['pythoncmd'] == sys.executable:
        if env['IS_WINDOWS']:
            python_inc_path=sysconfig.get_python_inc()
            python_lib_path=os.path.join(sysconfig.get_config_var('prefix'), 'libs')
            python_libs=['python%s%s'%(sys.version_info[0], sys.version_info[1])]
            verstring=".".join([str(i) for i in sys.version_info[:3]])
        else:
            (python_lib_path, python_libs,verstring, python_inc_path)=call_python_config()

    # if we want to use a python other than the one scons is running
    # Note: we assume scons is running python 2 in the following.
    else:
        (python_lib_path, python_libs,verstring, python_inc_path)=call_python_config(env['pythoncmd'])
    env['python_version'] = verstring
    ispython3 = (verstring[0] == '3')
    if ispython3:
        env.Append(CPPDEFINES=['ESPYTHON3'])
    env['buildvars']['python_version'] = verstring
    env['buildvars']['python'] = env['pythoncmd']
    # Check for an override from the config file.
    # Ideally, this should be automatic but we need to deal with the case
    # where python is not in its INSTALL directory
    if env['pythonlibpath'] != '':
        python_lib_path = env['pythonlibpath']

    if env['pythonincpath'] != '':
        python_inc_path = env['pythonincpath']

    if env['pythonlibname'] != '':
        python_libs = env['pythonlibname']

    conf = Configure(env.Clone())

    if env['sysheaderopt'] == '':
        conf.env.AppendUnique(CPPPATH = [python_inc_path])
    else:
        conf.env.Append(CCFLAGS = [env['sysheaderopt'], python_inc_path])

    conf.env.AppendUnique(LIBPATH = [python_lib_path])
    conf.env.AppendUnique(LIBS = python_libs)
    # The wrapper script needs to find the libs
    conf.env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], python_lib_path)

    if not conf.CheckCXXHeader('Python.h'):
        print("Cannot find python include files (tried 'Python.h' in directory %s)" % (python_inc_path))
        env.Exit(1)
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
    boost_inc_path,boost_lib_path=findLibWithHeader(env, env['boost_libs'], 'boost/python.hpp', env['boost_prefix'], lang='c++')
    if env['sysheaderopt'] == '':
        env.AppendUnique(CPPPATH = [boost_inc_path])
    else:
        # This is required because we can't -isystem /usr/include since it
        # breaks std includes
        if os.path.normpath(boost_inc_path) == '/usr/include':
            env.Append(CCFLAGS=[env['sysheaderopt'], os.path.join(boost_inc_path,'boost')])
        else:
            env.Append(CCFLAGS=[env['sysheaderopt'], boost_inc_path])

    env.AppendUnique(LIBPATH = [boost_lib_path])
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
    return env

def checkNumpy(env):
    if not detectModule(env, 'numpy'):
        print("Cannot import numpy. If it is installed try setting your PYTHONPATH and probably %s"%env['LD_LIBRARY_PATH_KEY'])
        env.Exit(1)

    ## check for numpy header (optional)
    conf = Configure(env.Clone())
    if conf.CheckCXXHeader(['Python.h','numpy/ndarrayobject.h']):
        conf.env.Append(CPPDEFINES = ['ESYS_HAVE_NUMPY_H'])
        conf.env['numpy_h']=True
    else:
        conf.env['numpy_h']=False

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
        env['warnings'].append("Cannot import scipy. NetCDF sources will not be available for inversions.")

    ######## pyproj
    if not detectModule(env, 'pyproj'):
        env['warnings'].append("Cannot import pyproj. Inversions may not work.")

    ######## gdal
    if not detectModule(env, 'gdal'):
        env['warnings'].append("Cannot import gdal. Inversions will not honour WKT coordinate system information.")

    ######## sympy
    if not detectModule(env, 'sympy'):
        env['warnings'].append("Cannot import sympy. Symbolic toolbox and nonlinear PDEs will not be available.")
    else:
        if env['pythoncmd'] is not None:
            get_external_python_sympy(env['pythoncmd'])
        else:
            import sympy as sp
            spVer=sp.__version__
            spl=spVer.split('.')
            if int(spl[0]) == 0 and int(spl[1]) < 7:
                env['sympy']=False
                env['warnings'].append("sympy version too old. Symbolic toolbox and nonlinear PDEs will not be available.")

    ######## gmshpy
    env['gmshpy'] = detectModule(env, 'gmshpy')

    return env

def checkForTrilinos(env):
    trilinos_inc_path=''
    trilinos_lib_path=''
    if env['trilinos']:
        havelibs = (len(env['trilinos_libs']) > 0)
        trilinos_inc_path,trilinos_lib_path=findLibWithHeader(env,
                env['trilinos_libs'], 'Tpetra_CrsMatrix.hpp',
                env['trilinos_prefix'], lang='c++', try_link=havelibs)
        if not havelibs:
            packages=['Tpetra','Kokkos','Belos','Amesos2','Ifpack2','MueLu']
            libs = []
            for pk in packages:
                # find out what libraries to link with...
                makefile = os.path.join(trilinos_inc_path, 'Makefile.export.%s'%pk)
                try:
                    for l in open(makefile, 'r').readlines():
                        if l.startswith("%s_LIBRARIES"%pk): # or l.startswith("Trilinos_TPL_LIBRARIES"):
                            lst = l.split('=')[1].strip().split()
                            lst = [e.replace('-l','',1) for e in lst]
                            libs += lst
                        elif l.startswith("%s_TPL_INCLUDE_DIRS"%pk):
                            lst = l.split('=')[1].strip().split()
                            lst = [e.replace('-I','',1) for e in lst]
                            env.AppendUnique(CPPPATH = lst)

                except Exception as e:
                    raise RuntimeError('Error reading Trilinos export Makefile\n%s'%(e))
            env['trilinos_libs'] = libs

        env.AppendUnique(CPPPATH = [trilinos_inc_path])
        env.AppendUnique(LIBPATH = [trilinos_lib_path])
        env.Append(CPPDEFINES = ['ESYS_HAVE_TRILINOS'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], trilinos_lib_path)
        env['buildvars']['trilinos_inc_path']=trilinos_inc_path
        env['buildvars']['trilinos_lib_path']=trilinos_lib_path
    env['buildvars']['trilinos']=int(env['trilinos'])
    return env

def checkOptionalLibraries(env):
    ######## netCDF
    netcdf_inc_path=''
    netcdf_lib_path=''
    if env['netcdf']:
        if env['netcdf']==4:
            env.Append(CPPDEFINES = ['NETCDF4'])
            netcdf_inc_path,netcdf_lib_path=findLibWithHeader(env, env['netcdf_libs'], 'ncVar.h', env['netcdf_prefix'], lang='c++')
        else:
            netcdf_inc_path,netcdf_lib_path=findLibWithHeader(env, env['netcdf_libs'], 'netcdfcpp.h', env['netcdf_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [netcdf_inc_path])
        env.AppendUnique(LIBPATH = [netcdf_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], netcdf_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_NETCDF'])
        env['buildvars']['netcdf_inc_path']=netcdf_inc_path
        env['buildvars']['netcdf_lib_path']=netcdf_lib_path
    env['buildvars']['netcdf']=int(env['netcdf'])

    ######## PAPI
    papi_inc_path=''
    papi_lib_path=''
    if env['papi']:
        papi_inc_path,papi_lib_path=findLibWithHeader(env, env['papi_libs'], 'papi.h', env['papi_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [papi_inc_path])
        env.AppendUnique(LIBPATH = [papi_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], papi_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_PAPI'])
        env['buildvars']['papi_inc_path']=papi_inc_path
        env['buildvars']['papi_lib_path']=papi_lib_path
    env['buildvars']['papi']=int(env['papi'])

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
                # try for clapack
                flavour = 'clapack'
                header = 'clapack.h'
                prefix = env['lapack_prefix']
                if len(env['lapack_libs']) == 0:
                    libs = ['lapack_atlas']
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
    if env['mpi']=='no':
        env['mpi']='none'

    env['usempi'] = env['mpi']!='none'
    mpi_inc_path=''
    mpi_lib_path=''
    if env['usempi']:
        mpi_inc_path,mpi_lib_path=findLibWithHeader(env, env['mpi_libs'], 'mpi.h', env['mpi_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [mpi_inc_path])
        env.AppendUnique(LIBPATH = [mpi_lib_path])
        env.AppendUnique(LIBS = env['mpi_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], mpi_lib_path)
        env.Append(CPPDEFINES = ['ESYS_MPI', 'MPI_NO_CPPBIND', 'MPICH_IGNORE_CXX_SEEK'])
        # NetCDF 4.1 defines MPI_Comm et al. if MPI_INCLUDED is not defined!
        # On the other hand MPT and OpenMPI don't define the latter so we have
        # to do that here
        if env['netcdf'] and env['mpi'] in ['MPT','OPENMPI']:
            env.Append(CPPDEFINES = ['MPI_INCLUDED'])

        if env['mpi'] == 'OPENMPI':
            # try to get version for correct launcher arguments
            try:
                p = Popen(['orterun', '-V'], stdout=PIPE, stderr=PIPE)
                o,e = p.communicate()
                try:
                    ver = e.split('\n')[0].split()[-1]
                except IndexError:
                    ver = o.split('\n')[0].split()[-1]
                if len(ver) > 0:
                    env['orte_version'] = ver
            except OSError:
                pass

        env['buildvars']['mpi_inc_path']=mpi_inc_path
        env['buildvars']['mpi_lib_path']=mpi_lib_path
    env['buildvars']['mpi']=env['mpi']

    ######## BOOMERAMG
    if env['mpi'] == 'none': env['boomeramg'] = False
    boomeramg_inc_path=''
    boomeramg_lib_path=''
    if env['boomeramg']:
        boomeramg_inc_path,boomeramg_lib_path=findLibWithHeader(env, env['boomeramg_libs'], 'HYPRE.h', env['boomeramg_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [boomeramg_inc_path])
        env.AppendUnique(LIBPATH = [boomeramg_lib_path])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], boomeramg_lib_path)
        env.Append(CPPDEFINES = ['ESYS_HAVE_BOOMERAMG'])
        env['buildvars']['boomeramg_inc_path']=boomeramg_inc_path
        env['buildvars']['boomeramg_lib_path']=boomeramg_lib_path
    env['buildvars']['boomeramg']=int(env['boomeramg'])

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

    ######## gmsh (for tests)
    env['gmsh'] = False
    if env['IS_WINDOWS']:
        try:
            p=Popen(['gmsh', '-info'], stderr=PIPE)
            _,e=p.communicate()
            if e.split().count("MPI"):
                env['gmsh']='m'
            else:
                env['gmsh']='s'
        except OSError:
            pass
    else:
        which = Popen(['which', 'gmsh'], stdout=PIPE)
        path,_ = which.communicate()
        if which.wait() == 0:
            cmd = ['ldd', path[:-1]]
            if env['IS_OSX']:
                cmd = ['otool','-L', path[:-1]]
            try:
                p=Popen(cmd, stdout=PIPE)
                gmshlibs,_ = p.communicate()
                env.Append(CPPDEFINES=['ESYS_HAVE_GMSH'])
                if p.returncode == 0 and 'libmpi' in gmshlibs:
                    env['gmsh'] = 'm'
                    env.Append(CPPDEFINES=['ESYS_GMSH_MPI'])
                else:
                    env['gmsh'] = 's'
            except OSError:
                pass
    
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
    return env

def checkPDFLatex(env):
    if 'PDF' in dir(env) and '.tex' in env.PDF.builder.src_suffixes(env):
        env['pdflatex']=True
    else:
        env['pdflatex']=False
    return env


