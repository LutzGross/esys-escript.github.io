# $Id$

"""

Tests importation of escript modules.

"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import os
import sys

from esys.escript import *
from esys.escript.linearPDEs import *
from esys.escript.pdetools import *

exit_code = 0

testdata = Data()
if testdata.isEmpty():
  print "Successfully created empty Data object"
else:
  exit_code = 1

sys.exit(exit_code)
