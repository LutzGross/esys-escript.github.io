
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""Classes and tools that form the basis of the escript system.
Specific solvers and domains are found in their respective packages."""


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os, sys

esys_packages = __name__.split('.')
esys_paths = os.path.dirname(os.path.abspath(__file__)).split(os.sep)
__buildinfo_file__ = None
if esys_paths[-(len(esys_packages)+1)] == 'site-packages':
    # check the library path for a pip install
    package_prefix = os.sep.join(esys_paths[:-len(esys_packages)])
    if os.name == 'nt':
        lib_env_name = 'PATH'
        cmd, env1, env2 = 'set', '%', '%'
    else:
        lib_env_name = 'LD_LIBRARY_PATH'
        cmd, env1, env2 = 'export', '$', ''
    esys_lib_path = os.sep.join([package_prefix, 'esys_escript_lib'])
    buildinfo_name = 'buildvars'
    __buildinfo_file__ = os.path.join(esys_lib_path, buildinfo_name)
    # check for the buildvars file
    if not os.path.exists(__buildinfo_file__):
        print('\ncreating {bv} file using {bv}.in template\n'.format(bv=buildinfo_name))
        with open(__buildinfo_file__+'.in', 'r') as f:
            lines = f.readlines()
        txt = ''.join(lines) % { 'python': sys.executable, 'prefix': package_prefix,
            'inc_prefix': '', 'lib_prefix': esys_lib_path }
        with open(__buildinfo_file__, 'w') as f:
            f.write(txt)
    lib_env_paths = os.environ.get(lib_env_name, '').split(os.pathsep)
    # re-start python if needed
    if esys_lib_path in lib_env_paths:
        print('\nusing escript library path: {}\n'.format(esys_lib_path))
    else:
        # set library path and restart python if not set
        new_lib_env_paths = os.sep.join([esys_lib_path] + lib_env_paths)
        print('\nwarning: escript library path not set')
        print('{} {}={}{}{}{}{}'.format(cmd, lib_env_name, esys_lib_path, os.pathsep, env1, lib_env_name, env2))
        env = os.environ
        env[lib_env_name] = os.pathsep.join(lib_env_paths+[esys_lib_path])
        if len(''.join(sys.argv)) > 0:
            print('re-starting python with new library path...\n')
            os.execve(sys.executable, [sys.executable]+sys.argv, env)
            sys.exit(0)
        else:
            # no args, don't restart an interactive session
            raise Exception('please exit python shell, set escript library path, then re-start')

from esys.escriptcore.escriptcpp import *
from esys.escriptcore.start import HAVE_SYMBOLS
from esys.escriptcore.util import *
from esys.escriptcore.nonlinearPDE import NonlinearPDE
from esys.escriptcore.datamanager import DataManager
if HAVE_SYMBOLS:
    from esys.escriptcore.symboliccore import *
from . import minimizer
from .domaincoupler import MPIDomainArray, DataCoupler
import logging
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)

__all__=[x for x in dir() if not x.startswith('internal_') and not x.startswith('Internal_') and not x.startswith('__') and not str(type(eval(x))).find('module')>=0]


