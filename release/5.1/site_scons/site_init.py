
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import sys, os, time, py_compile, re, subprocess
from SCons.Defaults import Chmod, Copy
from grouptest import *
from extractdebbuild import *

def findLibWithHeader(env, libs, header, paths, lang='c++', try_link=True):
    from SCons.Script.SConscript import Configure
    inc_path=''
    lib_path=''
    # 'paths' may be a prefix, so look for lib and include subdirectories
    if type(paths)==str:
        # find the header file first
        for i in 'include','include64','include32','inc':
            inc=os.path.join(paths, i)
            if os.path.isfile(os.path.join(inc, header)):
                inc_path=inc
                break
        if inc_path=='':
            raise RuntimeError('%s not found under %s'%(header,paths))

        # now try to find a lib directory
        for l in 'lib','lib64','lib32':
            lp=os.path.join(paths, l)
            if os.path.isdir(lp):
                lib_path=lp
                break
        if lib_path=='':
            raise RuntimeError('No lib directory found under %s'%paths)
    else:
        if os.path.isfile(os.path.join(paths[0], header)):
            inc_path=paths[0]
        else:
            raise RuntimeError('%s not found under %s'%(header,paths[0]))
        if os.path.isdir(paths[1]):
            lib_path=paths[1]
        else:
            raise RuntimeError('%s is not a valid path.'%paths[1])

    if try_link:
        # now try the library
        conf=Configure(env.Clone())
        conf.env.AppendUnique(CPPPATH = [inc_path])
        conf.env.AppendUnique(LIBPATH = [lib_path])
        if type(libs)==str: libs=[libs]
        if len(libs)==0: libs=['']
        # we can't check for each library by itself since they may depend on
        # each other, so we add all libraries to the link line and check only
        # for one
        conf.env.AppendUnique(LIBS = libs)
        if not conf.CheckLibWithHeader(libs[0], header, lang):
            conf.Finish()
            raise RuntimeError('Unable to link against %s (paths: %s, %s)'%(libs,inc_path,lib_path))

        conf.Finish()
    return inc_path, lib_path

def detectModule(env, module):
    from tempfile import TemporaryFile
    p=subprocess.call([env['pythoncmd'],'-c','import %s'%module], stderr=TemporaryFile())
    if p != 0:
        env[module] = False
        return False
    env[module] = True
    return True

def write_buildvars(env):
    buildvars=open(os.path.join(env['libinstall'], 'buildvars'), 'w')
    for k,v in sorted(env['buildvars'].items()):
        buildvars.write("%s=%s\n"%(k,v))
    buildvars.close()

def write_launcher(env):
    reps={'%n':'${ESCRIPT_NUM_NODES}', '%p':'${ESCRIPT_NUM_PROCS}',
          '%N':'${TOTPROC}', '%t':'${ESCRIPT_NUM_THREADS}', '%f':'${HOSTFILE}',
          '%h':'${HOSTLIST}', '%e':'${EXPORT_ENV}', '%b':'${EXEC_CMD}'}
    pre=env['prelaunch']
    cmd=env['launcher']
    post=env['postlaunch']
    # %b should be present in launcher at least
    if not '%b' in cmd:
        raise RuntimeError('option "launcher" must contain %b!')

    for k, v in reps.iteritems():
        pre = pre.replace(k, v)
        cmd = cmd.replace(k, v)
        post = post.replace(k, v)
    try:
        launchscript = os.path.join(env['bininstall'], 'run-escript')
        launcher=open(launchscript, 'w')
        if not env['stdlocationisprefix']:
            usestdlocation='0'
            stdlocation='/usr/lib/python-escript'
        else:
            usestdlocation='1'
            stdlocation=env['prefix']     
        for line in open('run-escript.in','r').readlines():
            s=line.replace('@@PRELAUNCH', pre).replace('@@LAUNCH', cmd).replace('@@POSTLAUNCH', post)
            s=s.replace('@@STDLOCATION', usestdlocation).replace('@@ESROOT',stdlocation)
            launcher.write(s)
        launcher.close()
        env.Execute(Chmod(launchscript, 0o755))
    except IOError:
        env['warnings'].append("Error attempting to write launcher script.")

def generateTestScripts(env, TestGroups):
    try:
        utest=open('utest.sh','w')
        utest.write(GroupTest.makeHeader(env['PLATFORM'], env['prefix'], False))
        for tests in TestGroups:
            utest.write(tests.makeString())
        utest.write(tests.makeFooter())
        utest.close()
        env.Execute(Chmod('utest.sh', 0o755))
        print("Generated utest.sh.")
        # This version contains only python tests - I want this to be usable
        # from a binary only install if you have the test files
        utest=open('itest.sh','w')
        utest.write(GroupTest.makeHeader(env['PLATFORM'], env['prefix'], True))
        for tests in TestGroups:
          if tests.exec_cmd=='$PYTHONRUNNER ':
            utest.write(tests.makeString())
        utest.write(tests.makeFooter())
        utest.close()
        env.Execute(Chmod('itest.sh', 0o755))
        print("Generated itest.sh.")        
    except IOError:
        env['warnings'].append("Error attempting to write unit test script(s).")

    # delete scripts upon cleanup
    env.Clean('target_init', 'utest.sh')
    env.Clean('target_init', 'itest.sh')

# Code to build .pyc from .py
def build_py(target, source, env):
    try:
       py_compile.compile(str(source[0]), str(target[0]), doraise=True)
       return 0
    except py_compile.PyCompileError, e:
       print(e)
       return 1


# Code to run unit_test executables
def runUnitTest(target, source, env):
  time_start = time.time()
  app = str(source[0].abspath)
  pn, sn= os.path.split(app)
  if not os.name== "nt":
     app = "cd "+pn+"; "+os.path.join(env['bininstall'], "run-escript")+" -bv "+os.path.join('.',sn)
  else:
      if env['usempi']:
          app = "cd %s & mpiexec -np %s -genvlist PYTHONPATH,OMP_NUM_THREADS,"\
            "FINLEY_TEST_DATA,PATH %s"\
            %(pn,env['ENV']['ESCRIPT_NUM_NODES'], sn)
      else:
           app = "cd "+ pn +" & "+sn
  print("Executing test: " + app)
  if not env.Execute(app):
    open(str(target[0]),'w').write("PASSED\n")
  else:
    return 1
  print("Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath))
  return None

def binpath(env, name=None):
    if not name:
        return env['bininstall']
    return os.path.join(env['bininstall'], name)

def runPyUnitTest(target, source, env): 
   time_start = time.time()
   app = str(source[0].abspath)
   pn, sn= os.path.split(app)
   if os.name=="nt":
       if env['usempi']:
           app = "cd %s & mpiexec -np %s -genvlist PYTHONPATH,OMP_NUM_THREADS,"\
              "FINLEY_TEST_DATA,PATH %s\pythonMPIredirect.exe %s"\
              %(pn,env['ENV']['ESCRIPT_NUM_NODES'],env['libinstall'],sn)
       else:
           app = "cd "+ pn +" & "+sys.executable + " " + sn
   else:
     skipfile = os.path.join(env['BUILD_DIR'], sn[:-3]) + ".skipped"
     failfile = os.path.join(env['BUILD_DIR'], sn[:-3]) + ".failed"
     try:
         os.unlink(skipfile)
     except Exception as e:
        pass
     app = "cd " + pn + "; " + binpath(env, "run-escript") + " -ov " + \
             str(env.Dir('#').abspath) + "/tools/testrunner.py -skipfile=" + \
             skipfile + " -failfile=" + failfile + " -exit " + sn
   print("Executing test: ",app)
   if env.Execute(app) == 0:
      open(str(target[0]),'w').write("PASSED\n")
   else:
     return 1
   print("Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath))
   return None

def runPyExample(target, source, env): 
   time_start = time.time()
   app = str(source[0].abspath)
   pn, sn= os.path.split(app)
   if os.name=="nt":
       if env['usempi']:
           app = "cd %s & mpiexec -np %s -genvlist PYTHONPATH,OMP_NUM_THREADS,"\
              "FINLEY_TEST_DATA,PATH %s\pythonMPIredirect.exe %s"\
              %(pn,env['ENV']['ESCRIPT_NUM_NODES'],env['libinstall'],sn)
       else:
           app = "cd "+ pn +" & "+sys.executable + " " + sn
   else:
    
     app = "cd "+pn+"; pwd; "+binpath(env, "run-escript")+" -ov "+sn
   print("Executing test: ",app)
   if env.Execute(app) == 0:
      open(str(target[0]),'w').write("PASSED\n")
   else:
     return 1
   print("Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath))
   return None

def eps2pdf(target, source, env):
#   if env.Execute("epstopdf "+str(source[0].abspath)+" -o "+str(target[0].abspath))!=0:
   if env.Execute("ps2pdf -dEPSCrop "+str(source[0].abspath)+" "+str(target[0].abspath))!=0:
       return 1
   return None

def effectiveName(inname):
    m=re.compile("^r1i[0-9]{1,2}n[0-9]{1,2}$")  # savanna names take the form r1i?n?
    if m.match(inname):
        return "savanna"
    return inname

def osxlib_dep_rewrite(libname, targetdir, env):
    if env.Execute("tools/libmover.sh %s %s"%(libname, targetdir)):
       return 1
    return None

def TristateVariable(key, help, default):
    """
    Modelled after SCons internal BoolVariable but allows three states
    (on=1, off=0, auto=-1)
    """
    on_strings = ('y', 'yes', 'true', 't', '1', 'on')
    off_strings = ('n', 'no', 'false', 'f', '0', 'off', 'none')
    auto_strings = ('a', 'auto', 'default', 'def', '-1', '')

    def _validator(key, val, env):
        if not env[key] in (1, 0, -1):
            raise SCons.Errors.UserError(
                    'Invalid value for tristate option %s: %s' % (key, env[key]))

    def _converter(val):
        lval = val.lower()
        if lval in on_strings: return 1
        if lval in off_strings: return 0
        if lval in auto_strings: return -1
        raise ValueError("Invalid value for tristate option: %s" % val)

    return (key, '%s (yes|no|auto)' % help, default, _validator, _converter)

