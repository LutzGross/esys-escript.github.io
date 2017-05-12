
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# This file includes setup and tweaks that are required for escript core packages
# No escript packages should be imported here


try:
    # this is required so newer intel MKL libraries find necessary symbols
    import ctypes, sys
    sys.setdlopenflags(sys.getdlopenflags()|ctypes.RTLD_GLOBAL)
except:
    pass

try:
    import sympy as sp
    spVer=sp.__version__
    spl=spVer.split('.')
    if int(spl[0]) == 0 and int(spl[1]) < 7:
        HAVE_SYMBOLS=False
    else:
        HAVE_SYMBOLS=True
except ImportError:
    HAVE_SYMBOLS=False

# To have this function called automatically
def escriptLogMemoryStatusNow(prefix='memescript'):
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
  


