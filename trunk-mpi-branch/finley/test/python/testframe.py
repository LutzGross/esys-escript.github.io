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
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle, JoinFaces, Brick

import numarray
FINLEY_TEST_MESH_PATH="data_meshes/"

NE=6 # number of element in each spatial direction (must be even)

class Test_X(unittest.TestCase):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.order=1
        d1 = Rectangle(n0=NE/2+1,n1=NE,l0=0.5,order=1)
        d2 = Rectangle(n0=NE/2,n1=NE,l0=0.5,order=1)
        d2.setX(d2.getX()+[0.5,0.])
        self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.order
        del self.domain
   def test_integrate_onFunctionOnBoundary_fromData_ReducedContinuousFunction_rank0(self):
      """
      tests integral of rank 0 Data on the FunctionOnBoundary

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports integral on FunctionOnBoundary
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnBoundary(self.domain)
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Data(0,(),w)
      if dim==2:
        arg=(-0.0177156089276)*x[0]+(-1.07750293477)*x[1]
        ref=(0.674554765151)*(1+2.*(dim-1.)/(o+1.))+(-1.76977330884)*dim
      else:
        arg=(0.304688056778)*x[0]+(0.548485298428)*x[1]+(0.672370309114)*x[2]
        ref=(0.0121419382123)*(1+2.*(dim-1.)/(o+1.))+(1.51340172611)*dim
      res=integrate(arg,where=w_ref)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")




if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_X))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


