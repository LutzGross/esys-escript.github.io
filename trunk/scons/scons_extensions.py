# Extensions to Scons

import py_compile
import sys
import os

# Code to build .pyc from .py
def build_py(target, source, env):
  py_compile.compile(str(source[0]), str(target[0]))
  return None

# Code to run unit_test executables
def runUnitTest(target, source, env):
  app = str(source[0].abspath)
  if not os.system(app):
    open(str(target[0]),'w').write("PASSED\n")
  else:
    return 1
  return None
