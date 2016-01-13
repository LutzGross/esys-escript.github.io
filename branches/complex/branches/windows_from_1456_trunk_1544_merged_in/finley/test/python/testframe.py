#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

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

NE=1 # number of element in each spatial direction (must be even)

class Test_X(unittest.TestCase):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   DEBUG=True
   def setUp(self):
        self.domain = Rectangle(n0=NE,n1=NE,l0=0.5,order=1)

   def test_setCoefficient_y_reduced_Scalar_using_y(self):
        d=self.domain.getDim()
        mypde=LinearPDE(self.domain,numSolutions=3,debug=self.DEBUG)
        mypde.setValue(y=Scalar(1.,ReducedFunctionOnBoundary(self.domain)))
        coeff=mypde.getCoefficientOfGeneralPDE("y_reduced")
        self.failUnlessEqual((coeff.getShape(),coeff.getFunctionSpace(),mypde.getNumEquations()),((),ReducedFunctionOnBoundary(self.domain),1))

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_X))
   s=unittest.TextTestRunner(verbosity=2).run(suite)


