# $Id$

import os
import sys

from escript.escript import *
from escript.linearPDEs import *
from escript.pdetools import *

"""

Tests importation of escript modules.

"""

exit_code = 0

testdata = Data()
if testdata.isEmpty():
  print "Successfully created empty Data object"
else:
  exit_code = 1

sys.exit(exit_code)
