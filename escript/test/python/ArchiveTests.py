
import sys
import unittest
import os

from esys.escript import *
from esys import bruce
from esys import finley

import numarray
from numarray import array,Float64,ones,greater

"""

Miscellaneous escript/Data tests.

Version $Id: MiscTests.py 153 2005-10-25 01:51:20Z jgs $

"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__licence__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licences/osl-3.0.php"""
#
# ==============================================================

print "\n\n"

mshList=(bruce.Rectangle(),
         bruce.Brick(),
         finley.Rectangle(2, 5, 1, l0 = 7.0, l1 = 11.0),
         finley.Brick(2, 5, 7, 1, l0 = 11.0, l1 = 13.0, l2 = 17.0),
         finley.Rectangle(2, 5, 2, l0 = 7.0, l1 = 11.0),
         finley.Brick(2, 5, 7, 2, l0 = 11.0, l1 = 13.0, l2 = 17.0))

for msh in mshList:

  print "\nX -- Continuous:"
  archDataX = msh.getX()
  archDataX.archiveData("data-archiveX")
  exDataX=Data()
  exDataX.extractData("data-archiveX",ContinuousFunction(msh))

  diff = archDataX - exDataX
  (infdiff, supdiff) = (inf(diff), sup(diff))
  if infdiff != 0 or supdiff != 0:
    print "*** ERROR: Data value discrepancies %f < X < %f" \
          % (infdiff, supdiff)
    sys.exit(1)

  exDataX.archiveData("data-archive2X");

  print "\nDataExpanded:"
  archDataE=Data([[1.00001],[2.00001]],Function(msh),True)
  archDataE.archiveData("data-archiveE")
  exDataE=Data()
  exDataE.extractData("data-archiveE",Function(msh))
  exDataE.archiveData("data-archive2E");

  print "\nDataTagged:"
  archDataT=Data([[1.00001],[2.00001]],Function(msh))
  archDataT.tag()
  archDataT.archiveData("data-archiveT")
  exDataT=Data()
  exDataT.extractData("data-archiveT",Function(msh))
  exDataT.archiveData("data-archive2T");

  print "\nDataConstant:"
  archDataC=Data([1.00001], Function(msh))
  archDataC.archiveData("data-archiveC")
  exDataC=Data()
  exDataC.extractData("data-archiveC",Function(msh))
  exDataC.archiveData("data-archive2C");

  print "\nDataEmpty:"
  archDataM=Data()
  archDataM.archiveData("data-archiveE")
  exDataM=Data()
  exDataM.extractData("data-archiveE",FunctionSpace())
  exDataM.archiveData("data-archive2E")

sys.exit(0)
# end
