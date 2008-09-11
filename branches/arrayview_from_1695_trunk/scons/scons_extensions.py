import sys, os, time, glob, fnmatch, types, py_compile

from SCons.Script.SConscript import SConsEnvironment

# Code to build .pyc from .py
def build_py(target, source, env):
    py_compile.compile(str(source[0]), str(target[0]))
    return 0

# Code to run unit_test executables
def runUnitTest(target, source, env):
  time_start = time.time()
  app = str(source[0].abspath)
  if env['usempi']: app = env['mpi_run'] + ' ' + app
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
   if env['usempi']:
     app = env['mpi_run'] +' lib/pythonMPI ' + app
   else:
     app = sys.executable + " " + app
   print "Executing test: " + app
   if env.Execute(app) == 0:
      open(str(target[0]),'w').write("PASSED\n")
   else:
     return 1
   print "Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath)
   return None

