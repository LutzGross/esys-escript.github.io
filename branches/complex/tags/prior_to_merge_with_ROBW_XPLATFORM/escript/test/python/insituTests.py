"""

Random interactive escript/Data tests.

Version $Id$

"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import os
import sys
import unittest
import math

from esys.escript import *
from esys import bruce

arglist = [ \
3.0, \
[3,4], \
[[1,2],[3,4]], \
[[15,8],[12,8]], \
[[[15,8],[12,8]],[[-9,9],[13,8]]] \
]

#
# ==============================================================

class escriptTestCase(unittest.TestCase):

  def setUp(self):
    self.msh=bruce.Rectangle()

  def prepareArg(self,val,ex,wh):
    if ex=="Expanded":
        exx=True
    else:
        exx=False
    out=Data(val,what=wh,expand=exx)
    if ex=="Tagged":
        out.tag()
    return out

  def testAlina1(self):
    P1 = 10.0
    assert log10(P1) == math.log10(10.0)
    assert log(P1) == math.log(10.0,math.e)

  def testAlina2(self):
    P = 10.0*Scalar(1.0, ContinuousFunction(self.msh))
    assert log10(P).convertToNumArray()[0] == math.log10(10.0)
    assert log(P).convertToNumArray()[0] == math.log(10.0,math.e)

  def testLog(self):
    for wh in [ContinuousFunction(self.msh),Function(self.msh)]:
      for ex in ["Constant","Tagged","Expanded"]:
        for a in arglist:
          #print "\n", ex, a, "==>"
          arg=self.prepareArg(a,ex,wh)
          #print "\nlog"
          result = arg._log10()

  def testLn(self):
    for wh in [ContinuousFunction(self.msh),Function(self.msh)]:
      for ex in ["Constant","Tagged","Expanded"]:
        for a in arglist:
          #print "\n", ex, a, "==>"
          arg=self.prepareArg(a,ex,wh)
          #print "\nln"
          result = arg._log()

  def testEmptyOp(self):
    emptyData=Data()
    emptyData._sin()
    emptyData2=Data()
    emptyData3=emptyData+emptyData2

if __name__ == '__main__':
  suite=unittest.TestSuite()
  suite.addTest(unittest.makeSuite(escriptTestCase))
  unittest.TextTestRunner(verbosity=2).run(suite)

sys.exit(0)
# end
