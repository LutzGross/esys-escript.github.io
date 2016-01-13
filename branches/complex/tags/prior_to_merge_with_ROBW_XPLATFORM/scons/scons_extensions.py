
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

# code to build epydoc docs
def build_epydoc(target, source, env):
    # get where I am currently, just as a reference
    pwd = os.getcwd()

    # get the full path of the runepydoc script
    runepydoc = str(source[0].abspath)

    # use this path to work out where the doc directory is
    dirs = runepydoc.split('/')
    dirs = dirs[:-3] # trim the last two entries: this is now the doc dir path
    docdir = '/'.join(dirs) # this is the backwards python way to do it
    # (I'm feeling in a perl mood today...)

    # change into the relevant dir
    os.chdir(docdir)

    # run the epydoc script
    if not os.system(runepydoc):
	os.chdir(pwd)
	open(str(target[0]), 'w').write("Documentation built\n")
    else:
	return 1
    return None

# build doxygen docs
def build_doxygen(target, source, env):
    # get where I am currently, just as a reference
    pwd = os.getcwd()

    # get the full path of the rundoxygen script
    rundoxygen = str(source[0].abspath)

    # use this path to work out where the doc directory is
    dirs = rundoxygen.split('/')
    dirs = dirs[:-2] # trim the last two entries: this is now the doc dir path
    docdir = '/'.join(dirs) # this is the backwards python way to do it
    # (I'm feeling in a perl mood today...)

    # change into the relevant dir
    os.chdir(docdir)

    # run the doxygen script
    if not os.system(rundoxygen):
	os.chdir(pwd)
	open(str(target[0]), 'w').write("Documentation built\n")
    else:
	return 1
    return None
