
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
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

def escriptOnExitProfiling():
    if os.name=='posix':
	pid=os.getpid()
	os.system("cat /proc/%d/status > memescript.%d"%(pid,pid))
	
if 'escriptExitProfiling' in os.environ:
	atexit.register(escriptOnExitProfiling)
