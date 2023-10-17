
"""Tool-specific initialization for clang-mp vi macport.
"""

# Based on SCons/Tool/gcc.py by Paweł Tomulik 2014 as a separate tool.

import os
import re
import subprocess

import SCons.Util
import SCons.Tool.cc
from SCons.Tool.clangCommon import get_clang_install_dirs
from SCons.Tool.MSCommon import msvc_setup_env_once

def generate(env, **kwargs):
    """Add Builders and construction variables for clang to an Environment."""
    compilers = ['clang++-mp-' + kwargs.get('version', '10')]
    SCons.Tool.cc.generate(env)
    if env['PLATFORM'] == 'win32':
        # Ensure that we have a proper path for clang
        clang = SCons.Tool.find_program_path(env, compilers[0], 
                                             default_paths=get_clang_install_dirs(env['PLATFORM']))
        if clang:
            clang_bin_dir = os.path.dirname(clang)
            env.AppendENVPath('PATH', clang_bin_dir)

            # Set-up ms tools paths
            msvc_setup_env_once(env)


    #print(env.Detect(compilers), compilers, env['CLANGMP_VERSION'])
    env['CXX'] = env.Detect(compilers) or compilers[0]
    if env['PLATFORM'] in ['cygwin', 'win32']:
        env['SHCXXFLAGS'] = SCons.Util.CLVar('$CXXFLAGS')
    else:
        env['SHCXXFLAGS'] = SCons.Util.CLVar('$CXXFLAGS -fPIC')

    # determine compiler version
    if env['CXX']:
        #pipe = SCons.Action._subproc(env, [env['CC'], '-dumpversion'],
        pipe = SCons.Action._subproc(env, [env['CXX'], '--version'],
                                     stdin='devnull',
                                     stderr='devnull',
                                     stdout=subprocess.PIPE)
        if pipe.wait() != 0: return
        # clang -dumpversion is of no use
        with pipe.stdout:
            line = pipe.stdout.readline()
        line = line.decode()
        match = re.search(r'clang +version +([0-9]+(?:\.[0-9]+)+)', line)
        if match:
            env['CXXVERSION'] = match.group(1)

def exists(env):
    return env.Detect(compilers)
