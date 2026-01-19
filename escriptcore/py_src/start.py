
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

"""
Startup configuration and environment setup for escript core packages.

This module performs initialization tasks required before other escript
packages are loaded, including setting library load flags for Intel MKL
compatibility and optional memory profiling setup.

:note: No escript packages should be imported in this module to prevent
       circular import issues.
"""


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

# This file includes setup and tweaks that are required for escript core packages
# No escript packages should be imported here

from esys.escriptcore.escriptcpp import hasFeature
try:
    # this is required so newer intel MKL libraries find necessary symbols
    import ctypes, sys
    sys.setdlopenflags(sys.getdlopenflags()|ctypes.RTLD_GLOBAL)
except:
    pass

if hasFeature('sympy'):

    try:
        import sympy as sp
        HAVE_SYMBOLS=True
    except ImportError:
        HAVE_SYMBOLS=False
else:
    HAVE_SYMBOLS=False

# To have this function called automatically
def escriptLogMemoryStatusNow(prefix='memescript'):
    """
    Logs the current memory status of the escript process to a file.

    On POSIX systems, this reads /proc/<pid>/status and writes it to a file
    named <prefix>.<pid>. This is useful for memory profiling and debugging.

    :param prefix: prefix for the output filename
    :type prefix: ``str``
    """
    import os
    if os.name=='posix':
        pid=os.getpid()
        os.system("cat /proc/%d/status > %s.%d"%(pid,prefix,pid))
    
try:
    import os
    if 'escriptExitProfiling' in os.environ:
        import atexit
        atexit.register(escriptLogMemoryStatusNow)
except:
    pass
  


