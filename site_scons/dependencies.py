
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

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import os, re, sys
from distutils import sysconfig
from subprocess import PIPE, Popen
from SCons.Script.SConscript import Configure
from site_init import findLibWithHeader, detectModule

def checkCompiler(env):
    conf = Configure(env.Clone())
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

    conf.env['buildvars']['cc']=conf.env['CC']
    conf.env['buildvars']['cxx']=conf.env['CXX']

    if conf.CheckFunc('gethostname'):
        conf.env.Append(CPPDEFINES = ['HAVE_GETHOSTNAME'])

    return conf.Finish()

def checkPython(env):
    # First we check to see if the config file has specified
    # where to find the file. Ideally, this should be automatic
    # but we need to deal with the case where python is not in its INSTALL
    # directory.
    # Use the python scons is running
    if env['pythoncmd']=='python':
        python_inc_path=sysconfig.get_python_inc()
        if env['IS_WINDOWS']:
            python_lib_path=os.path.join(sysconfig.get_config_var('prefix'), 'libs')
        elif env['PLATFORM']=='darwin':
            python_lib_path=sysconfig.get_config_var('LIBPL')
        else:
            python_lib_path=sysconfig.get_config_var('LIBDIR')

        #python_libs=[sysconfig.get_config_var('LDLIBRARY')] # only on linux
        if env['IS_WINDOWS']:
            python_libs=['python%s%s'%(sys.version_info[0], sys.version_info[1])]
        else:
            python_libs=['python'+sysconfig.get_python_version()]

        env['buildvars']['python']=sys.executable
        env['buildvars']['python_version']=str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2])

    #if we want to use a python other than the one scons is running
    else:
        initstring='from __future__ import print_function;from distutils import sysconfig;'
        if env['pythonlibname']!='':
            python_libs=env['pythonlibname']
        else:	# work it out by calling python    
            if ['IS_WINDOWS']:
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
        if ['IS_WINDOWS']:
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

        env['buildvars']['python']=env['pythoncmd']
        p=Popen([env['pythoncmd'], '-c', 'from __future__ import print_function;import sys;print(str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2]))'], stdout=PIPE)
        verstring=p.stdout.readline().strip()
        p.wait()
        env['buildvars']['python_version']=verstring

    # Check for an override from the config file.
    # Ideally, this should be automatic but we need to deal with the case
    # where python is not in its INSTALL directory
    if env['pythonlibpath']!='':
        python_lib_path=env['pythonlibpath']

    if env['pythonincpath']!='':
        python_inc_path=env['pythonincpath']

    conf = Configure(env.Clone())

    if env['sysheaderopt'] == '':
        conf.env.AppendUnique(CPPPATH = [python_inc_path])
    else:
        conf.env.Append(CCFLAGS = [env['sysheaderopt'], python_inc_path])

    conf.env.AppendUnique(LIBPATH = [python_lib_path])
    conf.env.AppendUnique(LIBS = python_libs)
    # The wrapper script needs to find the libs
    conf.env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], python_lib_path)

    if not conf.CheckCHeader('Python.h'):
        print("Cannot find python include files (tried 'Python.h' in directory %s)" % (python_inc_path))
        Exit(1)
    if not conf.CheckFunc('Py_Exit'):
        print("Cannot find python library method Py_Main (tried %s in directory %s)" % (python_libs, python_lib_path))
        Exit(1)

    return conf.Finish()

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
    try:
        for line in boosthpp:
            ver=re.match(r'#define BOOST_VERSION (\d+)',line)
            if ver:
                boostversion=ver.group(1)
    except StopIteration:
        pass
    boosthpp.close()
    env['buildvars']['boost_inc_path']=boost_inc_path
    env['buildvars']['boost_lib_path']=boost_lib_path
    env['buildvars']['boostversion']=boostversion
    return env

def checkNumpy(env):
    if not detectModule(env, 'numpy'):
        print("Cannot import numpy. If it is installed try setting your PYTHONPATH and probably %s"%env['LD_LIBRARY_PATH_KEY'])
        Exit(1)

    ## check for numpy header (optional)
    conf = Configure(env.Clone())
    if env['usepython3']:
        # FIXME: This is until we can work out how to make the checks in python 3
        conf.env['numpy_h']=False
    else:
        if conf.CheckCXXHeader(['Python.h','numpy/ndarrayobject.h']):
            conf.env.Append(CPPDEFINES = ['HAVE_NUMPY_H'])
            conf.env['numpy_h']=True
        else:
            conf.env['numpy_h']=False

    return conf.Finish()

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
    ######## sympy
    if not detectModule(env, 'sympy'):
        env['warnings'].append("Cannot import sympy. Symbolic toolbox and nonlinear PDEs will not be available.")

    ######## pyproj
    if not detectModule(env, 'pyproj'):
        env['warnings'].append("Cannot import pyproj. Inversions may not work.")

    ######## gdal
    if not detectModule(env, 'gdal'):
        env['warnings'].append("Cannot import gdal. Inversions will not honour WKT coordinate system information.")

    ######## scipy
    if not detectModule(env, 'scipy'):
        env['warnings'].append("Cannot import scipy. NetCDF sources will not be available for inversions.")

    return env

def checkOptionalLibraries(env):
    ######## netCDF
    netcdf_inc_path=''
    netcdf_lib_path=''
    if env['netcdf']:
        netcdf_inc_path,netcdf_lib_path=findLibWithHeader(env, env['netcdf_libs'], 'netcdf.h', env['netcdf_prefix'], lang='c++')
        env.AppendUnique(CPPPATH = [netcdf_inc_path])
        env.AppendUnique(LIBPATH = [netcdf_lib_path])
        env.AppendUnique(LIBS = env['netcdf_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], netcdf_lib_path)
        env.Append(CPPDEFINES = ['USE_NETCDF'])
        env['buildvars']['netcdf_inc_path']=netcdf_inc_path
        env['buildvars']['netcdf_lib_path']=netcdf_lib_path
    env['buildvars']['netcdf']=int(env['netcdf'])

    ######## PAPI
    papi_inc_path=''
    papi_lib_path=''
    if env['papi']:
        papi_inc_path,papi_lib_path=findLibWithHeader(env, env['papi_libs'], 'papi.h', env['papi_prefix'], lang='c')
        env.AppendUnique(CPPPATH = [papi_inc_path])
        env.AppendUnique(LIBPATH = [papi_lib_path])
        env.AppendUnique(LIBS = env['papi_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], papi_lib_path)
        env.Append(CPPDEFINES = ['BLOCKPAPI'])
        env['buildvars']['papi_inc_path']=papi_inc_path
        env['buildvars']['papi_lib_path']=papi_lib_path
    env['buildvars']['papi']=int(env['papi'])

    ######## MKL
    mkl_inc_path=''
    mkl_lib_path=''
    if env['mkl']:
        mkl_inc_path,mkl_lib_path=findLibWithHeader(env, env['mkl_libs'], 'mkl_solver.h', env['mkl_prefix'], lang='c')
        env.AppendUnique(CPPPATH = [mkl_inc_path])
        env.AppendUnique(LIBPATH = [mkl_lib_path])
        env.AppendUnique(LIBS = env['mkl_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], mkl_lib_path)
        env.Append(CPPDEFINES = ['MKL'])
        env['buildvars']['mkl_inc_path']=mkl_inc_path
        env['buildvars']['mkl_lib_path']=mkl_lib_path
    env['buildvars']['mkl']=int(env['mkl'])

    ######## UMFPACK
    umfpack_inc_path=''
    umfpack_lib_path=''
    if env['umfpack']:
        umfpack_inc_path,umfpack_lib_path=findLibWithHeader(env, env['umfpack_libs'], 'umfpack.h', env['umfpack_prefix'], lang='c')
        env.AppendUnique(CPPPATH = [umfpack_inc_path])
        env.AppendUnique(LIBPATH = [umfpack_lib_path])
        env.AppendUnique(LIBS = env['umfpack_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], umfpack_lib_path)
        env.Append(CPPDEFINES = ['UMFPACK'])
        env['buildvars']['umfpack_inc_path']=umfpack_inc_path
        env['buildvars']['umfpack_lib_path']=umfpack_lib_path
    env['buildvars']['umfpack']=int(env['umfpack'])

    ######## LAPACK
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
        env['buildvars']['lapack_inc_path']=lapack_inc_path
        env['buildvars']['lapack_lib_path']=lapack_lib_path
    env['buildvars']['lapack']=env['lapack']

    ######## Silo
    silo_inc_path=''
    silo_lib_path=''
    if env['silo']:
        silo_inc_path,silo_lib_path=findLibWithHeader(env, env['silo_libs'], 'silo.h', env['silo_prefix'], lang='c')
        env.AppendUnique(CPPPATH = [silo_inc_path])
        env.AppendUnique(LIBPATH = [silo_lib_path])
        # Note that we do not add the libs since they are only needed for the
        # weipa library and tools.
        #env.AppendUnique(LIBS = [env['silo_libs']])
        env['buildvars']['silo_inc_path']=silo_inc_path
        env['buildvars']['silo_lib_path']=silo_lib_path
    env['buildvars']['silo']=int(env['silo'])

    ######## VisIt
    visit_inc_path=''
    visit_lib_path=''
    if env['visit']:
        visit_inc_path,visit_lib_path=findLibWithHeader(env, env['visit_libs'], 'VisItControlInterface_V2.h', env['visit_prefix'], lang='c')
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
        mpi_inc_path,mpi_lib_path=findLibWithHeader(env, env['mpi_libs'], 'mpi.h', env['mpi_prefix'], lang='c')
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
        env['buildvars']['mpi_inc_path']=mpi_inc_path
        env['buildvars']['mpi_lib_path']=mpi_lib_path
    env['buildvars']['mpi']=env['mpi']

    ######## BOOMERAMG
    if env['mpi'] == 'none': env['boomeramg'] = False
    boomeramg_inc_path=''
    boomeramg_lib_path=''
    if env['boomeramg']:
        boomeramg_inc_path,boomeramg_lib_path=findLibWithHeader(env, env['boomeramg_libs'], 'HYPRE.h', env['boomeramg_prefix'], lang='c')
        env.AppendUnique(CPPPATH = [boomeramg_inc_path])
        env.AppendUnique(LIBPATH = [boomeramg_lib_path])
        env.AppendUnique(LIBS = env['boomeramg_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], boomeramg_lib_path)
        env.Append(CPPDEFINES = ['BOOMERAMG'])
        env['buildvars']['boomeramg_inc_path']=boomeramg_inc_path
        env['buildvars']['boomeramg_lib_path']=boomeramg_lib_path
    env['buildvars']['boomeramg']=int(env['boomeramg'])

    ######## ParMETIS
    if not env['usempi']: env['parmetis'] = False
    parmetis_inc_path=''
    parmetis_lib_path=''
    if env['parmetis']:
        parmetis_inc_path,parmetis_lib_path=findLibWithHeader(env, env['parmetis_libs'], 'parmetis.h', env['parmetis_prefix'], lang='c')
        env.AppendUnique(CPPPATH = [parmetis_inc_path])
        env.AppendUnique(LIBPATH = [parmetis_lib_path])
        env.AppendUnique(LIBS = env['parmetis_libs'])
        env.PrependENVPath(env['LD_LIBRARY_PATH_KEY'], parmetis_lib_path)
        env.Append(CPPDEFINES = ['USE_PARMETIS'])
        env['buildvars']['parmetis_inc_path']=parmetis_inc_path
        env['buildvars']['parmetis_lib_path']=parmetis_lib_path
    env['buildvars']['parmetis']=int(env['parmetis'])

    ######## gmsh (for tests)
    try:
        p=Popen(['gmsh', '-info'], stderr=PIPE)
        _,e=p.communicate()
        if e.split().count("MPI"):
            env['gmsh']='m'
        else:
            env['gmsh']='s'
    except OSError:
        env['gmsh']=False

    return env

def checkPDFLatex(env):
    if 'PDF' in dir(env) and '.tex' in env.PDF.builder.src_suffixes(env):
        env['pdflatex']=True
    else:
        env['pdflatex']=False
    return env


