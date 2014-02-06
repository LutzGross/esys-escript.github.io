
########################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for the linearPDE and pdetools on ripley

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import os
import unittest
from test_linearPDEs import Test_Poisson, Test_LinearPDE, Test_TransportPDE
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
from esys.ripley import Rectangle, Brick


try:
     RIPLEY_TEST_DATA=os.environ['RIPLEY_TEST_DATA']
except KeyError:
     RIPLEY_TEST_DATA='.'

NE=10 # number of element in each spatial direction (must be even)
mpiSize=getMPISizeWorld()

class Test_LinearPDEOnRipleyRect(Test_LinearPDE, Test_pdetools, Test_assemblage_2Do1, Test_TransportPDE):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

class Test_LinearPDEOnRipleyBrick(Test_LinearPDE, Test_pdetools, Test_assemblage_3Do1, Test_TransportPDE):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
            NX=x[0]
            NY=x[1]
            NZ=mpiSize//(x[0]*x[1])
            if NX*NY*NZ == mpiSize:
                break

        self.domain = Brick(n0=NE*NX-1, n1=NE*NY-1, n2=NE*NZ-1, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
        self.order = 1

    def tearDown(self):
        del self.domain

class Test_PoissonOnRipley(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
            NX=x
            NY=mpiSize//x
            if NX*NY == mpiSize:
                break
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    import sys
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_LinearPDEOnRipleyRect))
    suite.addTest(unittest.makeSuite(Test_LinearPDEOnRipleyBrick))
    suite.addTest(unittest.makeSuite(Test_PoissonOnRipley))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)

