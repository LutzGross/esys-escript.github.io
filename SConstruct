# top-level Scons configuration file for all esys13 modules
#
# set appropriate defaults for configuration variables
esysroot=str(Dir('.').abspath)
execfile(str(File(esysroot+"/scons/esys_options.py")))

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

#                  'doc/SConstruct']

SConscript(target_scripts, duplicate=0)
