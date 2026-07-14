import os, sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_linearPDEs import Test_Poisson, Test_LinearPDE_noLumping, Test_LameEquation, Test_Helmholtz
from esys.escript import getMPISizeWorld
from esys.oxley import Rectangle, Brick
NE=8
class Poisson2D(Test_Poisson):
    RES_TOL=1.e-7; ABS_TOL=1.e-8
    def setUp(self): self.domain=Rectangle(n0=NE, n1=NE, l0=1., l1=1.)
    def tearDown(self): del self.domain
class LinearPDE2D(Test_LinearPDE_noLumping, Test_LameEquation, Test_Helmholtz):
    RES_TOL=1.e-7; ABS_TOL=1.e-8
    def setUp(self): self.domain=Rectangle(n0=NE, n1=NE, l0=1., l1=1.); self.order=1
    def tearDown(self): del self.domain
class Poisson3D(Test_Poisson):
    RES_TOL=1.e-7; ABS_TOL=1.e-8
    def setUp(self): self.domain=Brick(n0=NE, n1=NE, n2=NE, l0=1., l1=1., l2=1.)
    def tearDown(self): del self.domain
class LinearPDE3D(Test_LinearPDE_noLumping, Test_LameEquation, Test_Helmholtz):
    RES_TOL=1.e-7; ABS_TOL=1.e-8
    def setUp(self): self.domain=Brick(n0=NE, n1=NE, n2=NE, l0=1., l1=1., l2=1.); self.order=1
    def tearDown(self): del self.domain
if __name__ == '__main__':
    suite = unittest.TestSuite()
    for c in [Poisson2D, LinearPDE2D, Poisson3D, LinearPDE3D]:
        suite.addTest(unittest.makeSuite(c))
    r = unittest.TextTestRunner(verbosity=1).run(suite)
    sys.exit(0 if r.wasSuccessful() else 1)
