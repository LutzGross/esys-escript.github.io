
#          Copyright 2006 by ACcESS MNRF                   
#                                                          
#              http://www.access.edu.au                    
#       Primary Business: Queensland, Australia            
#  Licensed under the Open Software License version 3.0    
#     http://www.opensource.org/licenses/osl-3.0.php       
#                                                          


# Extensions to Scons

import py_compile
import sys
import os
import time

# Code to build .pyc from .py
def build_py(target, source, env):
  py_compile.compile(str(source[0]), str(target[0]))
  return None

# Code to run unit_test executables
def runUnitTest(target, source, env):
  time_start = time.time()
  print "Executing test: " + str(source[0].abspath)
  app = str(source[0].abspath)
  if not env.Execute(app):
    open(str(target[0]),'w').write("PASSED\n")
  else:
    return 1
  print "Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath)
  return None

def runPyUnitTest(target, source, env): 
   time_start = time.time()
   print "Executing test: " + str(source[0].abspath)
   app = 'python '+'"'+str(source[0].abspath)+'"'
   if not env.Execute(app):
      open(str(target[0]),'w').write("PASSED\n")
   else:
     return 1
   print "Test execution time: ", round(time.time() - time_start, 1), " seconds wall time for " + str(source[0].abspath)
   return None
