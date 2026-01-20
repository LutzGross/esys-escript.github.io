
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""SCons.Tool.nvcc

Tool-specific initialization for NVIDIA CUDA Compiler.

There normally shouldn't be any need to import this module directly.
It will usually be imported through the generic SCons.Tool.Tool()
selection method.

This file copied with modifications from: http://www.scons.org/wiki/CudaTool
"""


import SCons.Tool
import SCons.Scanner.C
import SCons.Defaults
import subprocess
from tempfile import TemporaryFile

CUDASuffixes = ['.cu']

# make a CUDAScanner for finding #includes
# cuda uses the c preprocessor, so we can use the CScanner
CUDAScanner = SCons.Scanner.C.CScanner()

class ToolCudaWarning(SCons.Warnings.Warning):
    pass

class CudaCompilerNotFound(ToolCudaWarning):
    pass

SCons.Warnings.enableWarningClass(ToolCudaWarning)

def _detect(env):
    """ Try to detect the CUDA compiler """
    try:
        nvcc = env['NVCC']
    except KeyError:
        nvcc = env.WhereIs('nvcc')

    try:
        p=subprocess.call([nvcc,'-V'], stdout=TemporaryFile())
    except:
        raise SCons.Errors.StopError(CudaCompilerNotFound(
                "The NVIDIA CUDA compiler could not be found. Try setting "
                "'nvcc' in your options file or on the command line."))
    return nvcc

def add_common_nvcc_variables(env):
  """
  Add underlying common "NVIDIA CUDA compiler" variables that
  are used by multiple builders.
  """

  # "NVCC common command line"
  if not env.has_key('_NVCCCOMCOM'):
    # prepend -Xcompiler before each flag

    # these flags are common to both static and shared compilations
    env['_NVCCCOMCOM'] = '${_concat("-Xcompiler ", CPPFLAGS, "", __env__)} $_CPPDEFFLAGS $_CPPINCFLAGS'

    # wrap up all these environment variables inside -Xcompiler ""
    env['_NVCCWRAPCFLAGS'] =     '${_concat("-Xcompiler ", CFLAGS,     "", __env__)}'
    env['_NVCCWRAPSHCFLAGS'] =   '${_concat("-Xcompiler ", SHCFLAGS,   "", __env__)}'
    env['_NVCCWRAPCCFLAGS'] =    '${_concat("-Xcompiler ", CCFLAGS,    "", __env__)}'
    env['_NVCCWRAPSHCCFLAGS'] =  '${_concat("-Xcompiler ", SHCCFLAGS,  "", __env__)}'
    # XXX should these be wrapped as well?  not sure -jph
    #env['_NVCCWRAPCXXFLAGS'] =   '${_concat("-Xcompiler ", CXXFLAGS,   "", __env__)}'
    #env['_NVCCWRAPSHCXXFLAGS'] = '${_concat("-Xcompiler ", SHCXXFLAGS, "", __env__)}'

def generate(env):
  """
  Add Builders and construction variables for CUDA compilers to an Environment.
  """

  static_obj, shared_obj = SCons.Tool.createObjBuilders(env)

  for suffix in CUDASuffixes:
    # Add this suffix to the list of things buildable by Object
    static_obj.add_action('$CUDAFILESUFFIX', '$NVCCCOM')
    shared_obj.add_action('$CUDAFILESUFFIX', '$SHNVCCCOM')
    static_obj.add_emitter(suffix, SCons.Defaults.StaticObjectEmitter)
    shared_obj.add_emitter(suffix, SCons.Defaults.SharedObjectEmitter)

    # Add this suffix to the list of things scannable
    SCons.Tool.SourceFileScanner.add_scanner(suffix, CUDAScanner)

  add_common_nvcc_variables(env)

  # set the "CUDA Compiler Command" environment variable
  env['NVCC'] = _detect(env)
  env['SHNVCC'] = env['NVCC']
  
  # set the include path, and pass both c compiler flags and c++ compiler flags
  env['NVCCFLAGS'] = SCons.Util.CLVar('')
  env['SHNVCCFLAGS'] = SCons.Util.CLVar('') + ' -shared'
  
  # 'NVCC Command'
  env['NVCCCOM']   = '$NVCC -o $TARGET -c $NVCCFLAGS $_NVCCWRAPCFLAGS $_NVCCWRAPCCFLAGS $_NVCCCOMCOM $SOURCES'
  env['SHNVCCCOM'] = '$SHNVCC -o $TARGET -c $SHNVCCFLAGS $_NVCCWRAPSHCFLAGS $_NVCCWRAPSHCCFLAGS $_NVCCCOMCOM $SOURCES'
  
  # the suffix of CUDA source files is '.cu'
  env['CUDAFILESUFFIX'] = '.cu'


def exists(env):
  return _detect(env)

