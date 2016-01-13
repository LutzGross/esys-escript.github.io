# top-level Scons configuration file for all esys13 modules

#          Copyright 2006 by ACcESS MNRF                   
#                                                          
#              http://www.access.edu.au                    
#       Primary Business: Queensland, Australia            
#  Licensed under the Open Software License version 3.0    
#     http://www.opensource.org/licenses/osl-3.0.php       
#                                                          
#
#
# set appropriate defaults for configuration variables

esysroot=str(Dir('.').abspath)
execfile(str(File(esysroot+"/scons/esys_options.py")))

#
# call appropriate SConscripts

target_scripts = ['tools/CppUnitTest/SConstruct',
                  'esysUtils/SConstruct',
                  'escript/SConstruct',
                  'bruce/SConstruct',
                  'paso/SConstruct',
                  'finley/SConstruct',
                  'modellib/SConstruct',
                  'pyvisi/SConstruct']
                  # 'doc/SConstruct']
                  # 'doc/SConstruct']

SConscript(target_scripts, duplicate=0)
