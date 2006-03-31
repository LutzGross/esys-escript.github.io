# $Id:$
"""
frame to ran a single test out of the Test_util suite
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
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
       res=arg0+arg1
       ref=Data(numarray.array([5.4022518198315126, 4.2614777440630327]),self.functionspace)
       ref.setTaggedValue(1,numarray.array([3.7916899703627909, 2.650915894594311]))
       print arg0
       print arg1
       print res
       print ref
       self.failUnless(isinstance(res,Data),"wrong type of result.")
       self.failUnlessEqual(res.getShape(),(2,),"wrong shape of result.")
       self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result") 

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_mult_overloaded_constData_rank1_taggedData_rank0(self):
      arg0=Data(numarray.array([-1.8021292916234133, 0.52610779274754549]),self.functionspace)
      arg1=Data(-2.10992701792,self.functionspace)
      arg1.setTaggedValue(1,-0.103414543531)
      res=arg0*arg1
      ref=Data(numarray.array([3.8023612821737358, -1.1100490462541024]),self.functionspace)
      ref.setTaggedValue(1,numarray.array([0.18636637807672213, -0.054407197234984987]))
      self.failUnless(isinstance(res,Data),"wrong type of result.")
      self.failUnlessEqual(res.getShape(),(2,),"wrong shape of result.")
      self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_pow_overloaded_taggedData_rank1_taggedData_rank0(self):
      arg0=Data(numarray.array([2.3271463200997102, 1.1183883755594386]),self.functionspace)
      arg0.setTaggedValue(1,numarray.array([3.8629630953883574, 3.7380725394052305]))
      arg1=Data(4.96280973024,self.functionspace)
      arg1.setTaggedValue(1,1.353354768)
      res=arg0**arg1
      ref=Data(numarray.array([66.141826047243441, 1.7424328429281473]),self.functionspace)
      ref.setTaggedValue(1,numarray.array([6.2274713351642141, 5.9565601105055705]))
      self.failUnless(isinstance(res,Data),"wrong type of result.")
      self.failUnlessEqual(res.getShape(),(2,),"wrong shape of result.")
      self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
   def test_quotient_overloaded_constData_rank1_taggedData_rank0(self):
      arg0=Data(numarray.array([3.2585370968848757, 1.1454175877624291]),self.functionspace)
      arg1=Data(-2.93544378089,self.functionspace)
      arg1.setTaggedValue(1,-1.7765826732)
      res=arg0/arg1
      ref=Data(numarray.array([-1.1100662591795942, -0.39020252924586424]),self.functionspace)
      ref.setTaggedValue(1,numarray.array([-1.8341601244001815, -0.64473081103446983]))
      self.failUnless(isinstance(res,Data),"wrong type of result.")
      self.failUnlessEqual(res.getShape(),(2,),"wrong shape of result.")
      self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_util2))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
