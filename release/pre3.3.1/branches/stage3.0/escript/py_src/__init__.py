
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from escript import *
from util import *

import os
import atexit

# To have this function called automatically
def escriptLogMemoryStatusNow(prefix='memescript'):
    if os.name=='posix':
	pid=os.getpid()
	os.system("cat /proc/%d/status > %s.%d"%(pid,prefix,pid))
	
if 'escriptExitProfiling' in os.environ:
	atexit.register(escriptLogMemoryStatusNow)
