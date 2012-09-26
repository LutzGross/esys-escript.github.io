
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

try:
    # this is required so newer intel MKL libraries find necessary symbols
    import ctypes, sys
    sys.setdlopenflags(sys.getdlopenflags()|ctypes.RTLD_GLOBAL)
except:
    pass

try:
    import sympy
    HAVE_SYMBOLS=True
except ImportError:
    HAVE_SYMBOLS=False

#from .escript import *
from .escriptcpp import *
from .util import *
from .nonlinearPDE import NonlinearPDE
from .datamanager import DataManager

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

