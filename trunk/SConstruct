# top-level Scons configuration file for all esys13 modules

import os

#
# ensure correct versions of python and scons

EnsurePythonVersion(2,3)
EnsureSConsVersion(0,96)

#
# retreive command-line arguments if any
# these are passed through to SConscripts

if ARGUMENTS.get('libinstall',0):
  libinstall = ARGUMENTS.get('libinstall',0)
else:
  libinstall = Dir('#lib')
Export(["libinstall"])

if ARGUMENTS.get('options',0):
  options = ARGUMENTS.get('options',0)
else:
  options = None
Export(["options"])

if ARGUMENTS.get('debug',0):
  dodebug = 1
else:
  dodebug = 0
Export(["dodebug"])

if ARGUMENTS.get('usegcc',0):
  usegcc = 1
else:
  usegcc = 0
Export(["usegcc"])

#
# set and export esysroot

esysroot = Dir('#.')
Export(["esysroot"])

#
# call appropriate SConscripts

target_scripts = ['tools/CppUnitTest/SConstruct',
                  'tools/mmio/SConstruct',
                  'esysUtils/SConstruct',
                  'escript/SConstruct',
                  'bruce/SConstruct',
                  'paso/SConstruct',
                  'finley/SConstruct',
                  'modellib/SConstruct']

SConscript(target_scripts, duplicate=0)
