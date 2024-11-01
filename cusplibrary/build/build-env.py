
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################
from __future__ import print_function, division

EnsureSConsVersion(1,2)

import os

import inspect
import platform

def get_cuda_paths():
  """Determines CUDA {bin,lib,include} paths
  
  returns (bin_path,lib_path,inc_path)
  """

  # determine defaults
  if os.name == 'nt':
    bin_path = 'C:/CUDA/bin'
    lib_path = 'C:/CUDA/lib'
    inc_path = 'C:/CUDA/include'
  elif os.name == 'posix':
    bin_path = '/usr/bin'
    lib_path = '/usr/lib'
    inc_path = '/usr/include'
  else:
    raise ValueError, 'Error: unknown OS.  Where is nvcc installed?'
   
  if platform.machine()[-2:] == '64':
    lib_path += '64'

  # override with environement variables
  if 'CUDA_BIN_PATH' in os.environ:
    bin_path = os.path.abspath(os.environ['CUDA_BIN_PATH'])
  if 'CUDA_LIB_PATH' in os.environ:
    lib_path = os.path.abspath(os.environ['CUDA_LIB_PATH'])
  if 'CUDA_INC_PATH' in os.environ:
    inc_path = os.path.abspath(os.environ['CUDA_INC_PATH'])

  return (bin_path,lib_path,inc_path)

def get_mkl_paths():
  """Determines MKL {lib,include} paths
  
  returns (lib_path,inc_path)
  """

  arch64 = False
  if platform.machine()[-2:] == '64':
    arch64 = True

  if 'MKLROOT' not in os.environ:
    raise ValueError, "MKLROOT is not an environment variable"

  # determine defaults
  if os.name == 'nt':
    raise ValueError, "Intel MKL support for Windows not implemented."
  elif os.name == 'posix':
    lib_base = os.environ['MKLROOT'] + '/lib'
    dirs = os.listdir(lib_base)
    for dir in dirs :
        # select 64/32 bit MKL library path based on architecture
        if arch64 == True and dir.find('64') > -1 :
          lib_path = lib_base + '/' + dir
          break
        elif arch64 == False and dir.find('64') == -1 :
          lib_path = lib_base + '/' + dir
          break

    if lib_path == lib_base :
      raise ValueError, 'Could not find MKL library directory which matches the arctitecture.'

    inc_path = os.environ['MKLROOT'] + '/include'
  else:
    raise ValueError, 'Error: unknown OS.  Where is nvcc installed?'
   
  return (lib_path,inc_path)

def getTools():
  result = []
  if os.name == 'nt':
    result = ['default', 'msvc']
  elif os.name == 'posix':
    result = ['default', 'gcc']
  else:
    result = ['default']
  return result;


OldEnvironment = Environment;


# this dictionary maps the name of a compiler program to a dictionary mapping the name of
# a compiler switch of interest to the specific switch implementing the feature
gCompilerOptions = {
    'gcc' : {'warn_all' : '-Wall', 'warn_errors' : '-Werror', 'optimization' : '-O2', 'debug' : '-g',  'exception_handling' : '',      'omp' : '-fopenmp'},
    'g++' : {'warn_all' : '-Wall', 'warn_errors' : '-Werror', 'optimization' : '-O2', 'debug' : '-g',  'exception_handling' : '',      'omp' : '-fopenmp'},
    'cl'  : {'warn_all' : '/Wall', 'warn_errors' : '/WX',     'optimization' : '/Ox', 'debug' : ['/Zi', '-D_DEBUG', '/MTd'], 'exception_handling' : '/EHsc', 'omp' : '/openmp'}
  }


# this dictionary maps the name of a linker program to a dictionary mapping the name of
# a linker switch of interest to the specific switch implementing the feature
gLinkerOptions = {
    'gcc' : {'debug' : ''},
    'g++' : {'debug' : ''},
    'link'  : {'debug' : '/debug'}
  }


def getCFLAGS(mode, backend, warn, warnings_as_errors, hostspblas, b40c, CC):
  result = []
  if mode == 'release':
    # turn on optimization
    result.append(gCompilerOptions[CC]['optimization'])
  elif mode == 'debug':
    # turn on debug mode
    result.append(gCompilerOptions[CC]['debug'])
    result.append('-DTHRUST_DEBUG')
  # force 32b code on darwin
  if platform.platform()[:6] == 'Darwin':
    result.append('-m32')
  # build with B40C enabled 
  if b40c == True :
    result.append('-D__CUSP_USE_B40C__') 
  
  if CC == 'cl':
    result.append('/bigobj')

  # generate omp code
  if backend == 'omp':
    result.append(gCompilerOptions[CC]['omp'])

  if warn:
    # turn on all warnings
    result.append(gCompilerOptions[CC]['warn_all'])

  if warnings_as_errors:
    # treat warnings as errors
    result.append(gCompilerOptions[CC]['warn_errors'])

  # generate hostspblas code
  if hostspblas == 'mkl':
    result.append('-DINTEL_MKL_SPBLAS')

  return result


def getCXXFLAGS(mode, backend, warn, warnings_as_errors, hostspblas, b40c, CXX):
  result = []
  if mode == 'release':
    # turn on optimization
    result.append(gCompilerOptions[CXX]['optimization'])
  elif mode == 'debug':
    # turn on debug mode
    result.append(gCompilerOptions[CXX]['debug'])
  # enable exception handling
  result.append(gCompilerOptions[CXX]['exception_handling'])
  # force 32b code on darwin
  if platform.platform()[:6] == 'Darwin':
    result.append('-m32')
  # build with B40C enabled 
  if b40c is not None :
    result.append('-D__CUSP_USE_B40C__') 

  # generate omp code
  if backend == 'omp':
    result.append(gCompilerOptions[CXX]['omp'])

  if warn:
    # turn on all warnings
    result.append(gCompilerOptions[CXX]['warn_all'])

  if warnings_as_errors:
    # treat warnings as errors
    result.append(gCompilerOptions[CXX]['warn_errors'])

  # generate hostspblas code
  if hostspblas == 'mkl':
    result.append('-DINTEL_MKL_SPBLAS')

  return result


def getNVCCFLAGS(mode, backend, arch):
  result = ['-arch=' + arch]
  if mode == 'debug':
    # turn on debug mode
    # XXX make this work when we've debugged nvcc -G
    #result.append('-G')
    pass
  return result


def getLINKFLAGS(mode, backend, hostspblas, LINK):
  result = []
  if mode == 'debug':
    # turn on debug mode
    result.append(gLinkerOptions[LINK]['debug'])
  # force 32b code on darwin
  if platform.platform()[:6] == 'Darwin':
    result.append('-m32')

  # XXX make this portable
  if backend == 'ocelot':
    result.append(os.popen('OcelotConfig -l').read().split())

  if hostspblas == 'mkl':
    result.append('-fopenmp')

  return result


def Environment():
  # allow the user discretion to choose the MSVC version
  vars = Variables()
  if os.name == 'nt':
    vars.Add(EnumVariable('MSVC_VERSION', 'MS Visual C++ version', None, allowed_values=('8.0', '9.0', '10.0')))

  # add a variable to handle the device backend
  backend_variable = EnumVariable('backend', 'The parallel device backend to target', 'cuda',
                                  allowed_values = ('cuda', 'omp', 'ocelot'))
  vars.Add(backend_variable)

  # add a variable to handle RELEASE/DEBUG mode
  vars.Add(EnumVariable('mode', 'Release versus debug mode', 'release',
                        allowed_values = ('release', 'debug')))

  # add a variable to handle compute capability
  vars.Add(EnumVariable('arch', 'Compute capability code generation', 'sm_30',
                        allowed_values = ('sm_10', 'sm_11', 'sm_12', 'sm_13', 'sm_20', 'sm_21', 'sm_30', 'sm_35')))

  # add a variable to handle warnings
  if os.name == 'posix':
    vars.Add(BoolVariable('Wall', 'Enable all compilation warnings', 1))
  else:
    vars.Add(BoolVariable('Wall', 'Enable all compilation warnings', 0))

  # add a variable to treat warnings as errors
  vars.Add(BoolVariable('Werror', 'Treat warnings as errors', 0))

  # add a variable to filter source files by a regex
  vars.Add('tests', help='Filter test files using a regex')

  # add a variable to handle the device backend
  hostspblas_variable = EnumVariable('hostspblas', 'Host sparse math library', 'cusp',
                                  allowed_values = ('cusp', 'mkl'))
  vars.Add(hostspblas_variable)

  # add a variable to enable B40C support
  vars.Add(BoolVariable('b40c', 'Enable support for B40C', 0))

  # create an Environment
  env = OldEnvironment(tools = getTools(), variables = vars)

  # get the absolute path to the directory containing
  # this source file
  thisFile = inspect.getabsfile(Environment)
  thisDir = os.path.dirname(thisFile)

  # enable nvcc
  env.Tool('nvcc', toolpath = [os.path.join(thisDir)])

  # get the preprocessor define to use for the backend
  backend_define = { 'cuda' : 'THRUST_DEVICE_SYSTEM_CUDA', 'omp' : 'THRUST_DEVICE_SYSTEM_OMP', 'ocelot' : 'THRUST_DEVICE_SYSTEM_CUDA' }[env['backend']] 
  env.Append(CFLAGS = ['-DTHRUST_DEVICE_SYSTEM=%s' % backend_define])
  env.Append(CXXFLAGS = ['-DTHRUST_DEVICE_SYSTEM=%s' % backend_define])

  # get C compiler switches
  env.Append(CFLAGS = getCFLAGS(env['mode'], env['backend'], env['Wall'], env['Werror'], env['hostspblas'], env['b40c'], env.subst('$CC')))

  # get CXX compiler switches
  env.Append(CXXFLAGS = getCXXFLAGS(env['mode'], env['backend'], env['Wall'], env['Werror'], env['hostspblas'], env['b40c'], env.subst('$CXX')))

  # get NVCC compiler switches
  env.Append(NVCCFLAGS = getNVCCFLAGS(env['mode'], env['backend'], env['arch']))

  # get linker switches
  env.Append(LINKFLAGS = getLINKFLAGS(env['mode'], env['backend'], env['hostspblas'], env.subst('$LINK')))
   
  # silence unknown pragma warnings
  env.AppendUnique(CFLAGS = ['-Wno-unknown-pragmas','-Wno-unused-local-typedefs'])
  env.AppendUnique(CXXFLAGS = ['-Wno-unknown-pragmas','-Wno-unused-local-typedefs'])

  # get CUDA paths
  (cuda_exe_path,cuda_lib_path,cuda_inc_path) = get_cuda_paths()
  env.AppendUnique(LIBPATH = [cuda_lib_path])
  env.AppendUnique(CPPPATH = [cuda_inc_path])

  # link against backend-specific runtimes
  # XXX we shouldn't have to link against cudart unless we're using the
  #     cuda runtime, but cudafe inserts some dependencies when compiling .cu files
  # XXX ideally this gets handled in nvcc.py if possible
  env.AppendUnique(LIBS = ['cudart','stdc++','m'])

  if env['backend'] == 'ocelot':
    if os.name == 'posix':
      env.AppendUnique(LIBPATH = ['/usr/local/lib'])
    else:
      raise ValueError, "Unknown OS.  What is the Ocelot library path?"
  elif env['backend'] == 'omp':
    if os.name == 'posix':
      env.Append(LIBS = ['gomp'])
    elif os.name == 'nt':
      env.Append(LIBS = ['VCOMP'])
    else:
      raise ValueError, "Unknown OS.  What is the name of the OpenMP library?"

  if env['hostspblas'] == 'mkl':
    intel_lib = 'mkl_intel'
    if platform.machine()[-2:] == '64':
        intel_lib += '_lp64'
    
    (mkl_lib_path,mkl_inc_path) = get_mkl_paths()
    env.AppendUnique(CPPPATH = [mkl_inc_path])
    env.AppendUnique(LIBPATH = [mkl_lib_path])
    env.AppendUnique(LIBS = ['mkl_core', 'mkl_gnu_thread', intel_lib])

  # set thrust include path
  # this needs to come before the CUDA include path appended above,
  # which may include a different version of thrust
  env.Prepend(CPPPATH = os.path.dirname(thisDir))


  if 'THRUST_PATH' in os.environ:
      env.Prepend(CPPPATH = [os.path.abspath(os.environ['THRUST_PATH'])])

  # import the LD_LIBRARY_PATH so we can run commands which depend
  # on shared libraries
  # XXX we should probably just copy the entire environment
  if os.name == 'posix':
    if ('DYLD_LIBRARY_PATH' in os.environ) and (env['PLATFORM'] == "darwin") :
        env['ENV']['DYLD_LIBRARY_PATH'] = os.environ['DYLD_LIBRARY_PATH']
    elif 'LD_LIBRARY_PATH' in os.environ:
      env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']

  # generate help text
  Help(vars.GenerateHelpText(env))

  return env

