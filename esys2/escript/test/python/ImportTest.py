# $Id$

import os
import sys

#esys_root=os.getenv('ESYS_ROOT')
#sys.path.append(esys_root+'/escript/lib')
#sys.path.append(esys_root+'/escript/py_src')
                                                                                                                 
from escript.escript import *
from escript.linearPDEs import *
from escript.pdetools import *

"""

Tests importation of escript modules.

"""

testdata = Data()

if testdata.isEmpty():
  print "Successfully created empty Data object"

print "If an 8 is printed here => ", OPENINVENTOR, " <= util was imported OK."
