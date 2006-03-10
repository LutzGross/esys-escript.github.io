# $Id:$
"""
frame to ran a single test out of the Test_util suite
"""

import unittest
from esys.escript import *
from esys.finley import Rectangle
import numarray

class Test_util2(unittest.TestCase):
   RES_TOL=1.e-7
   def setUp(self):
       self.__dom =Rectangle(10,10,2)
       self.functionspace = FunctionOnBoundary(self.__dom) # due to a bug in escript python needs to hold a reference to the domain
   def test_add_overloaded_constData_rank1_taggedData_rank0(self):
       arg0=Data(numarray.array([4.5897569702707663, 3.4489828945022865]),self.functionspace)
       arg1=Data(0.812494849561,self.functionspace)
       arg1.setTaggedValue(1,-0.798066999908)
       print arg0
       print arg1
       res=arg0+arg1
       print res
       ref=Data(numarray.array([5.4022518198315126, 4.2614777440630327]),self.functionspace)
       ref.setTaggedValue(1,numarray.array([3.7916899703627909, 2.650915894594311]))
       self.failUnless(isinstance(res,Data),"wrong type of result.")
       self.failUnlessEqual(res.getShape(),(2,),"wrong shape of result.")
       self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result") 
if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_util2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
