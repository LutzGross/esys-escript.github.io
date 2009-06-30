
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import sys, os, time, glob, fnmatch, types, py_compile, re

from SCons.Script.SConscript import SConsEnvironment

# Code to build .pyc from .py
def build_py(target, source, env):
    py_compile.compile(str(source[0]), str(target[0]))
    return 0

# Code to run unit_test executables
def runUnitTest(target, source, env):
  time_start = time.time()
  app = str(source[0].abspath)
  pn, sn= os.path.split(app)
  if not os.name== "nt":
     app = "cd "+pn+"; "+os.path.join(env['bininstall'],"escript")+" -bv "+sn
  else:
      if env['usempi']:
          app = "cd %s | mpiexec -np %s -genvlist PYTHONPATH,OMP_NUM_THREADS,"\
            "FINLEY_TEST_DATA,PYVISI_TEST_DATA_ROOT,PYVISI_WORKDIR,PATH %s"\
            %(pn,env['ENV']['ESCRIPT_NUM_NODES'], sn)
      else:
           app = "cd "+ pn +" | "+ sn
  print "Executing test: " + app
  if not env.Execute(app):
    open(str(target[0]),'w').write("PASSED\n")
  else:
    return 1
  print "Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath)
  return None

def runPyUnitTest(target, source, env): 
   time_start = time.time()
   app = str(source[0].abspath)
   pn, sn= os.path.split(app)
   if os.name== "nt":
       if env['usempi']:
           app = "cd %s | mpiexec -np %s -genvlist PYTHONPATH,OMP_NUM_THREADS,"\
              "FINLEY_TEST_DATA,PYVISI_TEST_DATA_ROOT,PYVISI_WORKDIR,PATH %s\pythonMPIredirect.exe %s"\
              %(pn,env['ENV']['ESCRIPT_NUM_NODES'],env['libinstall'],sn)
       else:
           app = "cd "+ pn +" | "+sys.executable + " " + sn
   else:
     app = "cd "+pn+"; "+os.path.join(env['bininstall'],"escript")+" -ov "+sn
   print "Executing test: ",app
   if env.Execute(app) == 0:
      open(str(target[0]),'w').write("PASSED\n")
   else:
     return 1
   print "Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath)
   return None

def eps2pdf(target, source, env):
   if env.Execute("epstopdf "+str(source[0].abspath)+" -o "+str(target[0].abspath))!=0:
	   return 1
   return None

def effectiveName(inname):
   m=re.compile("^r1i[0-9]{1,2}n[0-9]{1,2}$")	# savanna names take the form r1i?n?
   if m.match(inname):
	return "service0"
   return inname
