# Extensions to Scons

import py_compile
import sys

# Code to build .pyc from .py
def build_py(target, source, env):
  py_compile.compile(str(source[0]), str(target[0]))
  return None
