"""

Random interactive escript/Data tests.

Version $Id$

"""

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
    return out

  def testAlina1(self):
    P1 = 10.0
    assert log(P1) == math.log10(10.0)
    assert ln(P1) == math.log(10.0,math.e)

  def testAlina2(self):
    P = 10.0*Scalar(1.0, ContinuousFunction(self.msh))
    assert log(P).convertToNumArray()[0] == math.log10(10.0)
    assert ln(P).convertToNumArray()[0] == math.log(10.0,math.e)

  def testLog(self):
    for wh in [ContinuousFunction(self.msh),Function(self.msh)]:
      for ex in ["Constant","Expanded"]:
        for a in arglist:
          #print "\n", ex, a, "==>"
          arg=self.prepareArg(a,ex,wh)
          #print "\nlog"
          result = arg.log()

  def testLn(self):
    for wh in [ContinuousFunction(self.msh),Function(self.msh)]:
      for ex in ["Constant","Expanded"]:
        for a in arglist:
          #print "\n", ex, a, "==>"
          arg=self.prepareArg(a,ex,wh)
          #print "\nln"
          result = arg.ln()

if __name__ == '__main__':
  suite=unittest.TestSuite()
  suite.addTest(unittest.makeSuite(escriptTestCase))
  unittest.TextTestRunner(verbosity=2).run(suite)

sys.exit(0)
# end
