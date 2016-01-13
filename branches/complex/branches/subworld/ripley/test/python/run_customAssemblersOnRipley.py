
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.ripley import Rectangle, Brick
from esys.escript.linearPDEs import LameEquation

EXPANDED, SCALAR, CONSTANT = range(3)

class RipleyCustomAssemblerTestBase(unittest.TestCase):
    def run_lame(self, fast, test_type, mu=3, lamb=50):
        d=self.domain.getDim()
        mypde=LameEquation(self.domain, useFast=fast)
        cf=ContinuousFunction(self.domain)
        x=cf.getX()
        u_ex=x
        msk=Vector(0.,cf)
        for i in range(d):
            msk[i]=whereZero(x[i])
        if test_type != CONSTANT:
            mu = Scalar(mu, cf)
            lamb = Scalar(lamb, cf)
            if test_type == EXPANDED:
                mu.expand()
                lamb.expand()
        mypde.setValue(q=msk,r=u_ex,lame_mu=mu,lame_lambda=lamb,f=(2*3+50*d)*FunctionOnBoundary(self.domain).getNormal())
        return mypde.getSolution()

    def test_lameExpanded(self):
        #check default and lame assemblers agree
        default = self.run_lame(False, EXPANDED)
        lame = self.run_lame(True, EXPANDED)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "expanded data"%self.domain.getDim())
        #reverse order, ensure default assembler still operational
        lame = self.run_lame(True, EXPANDED, mu=7, lamb=40)
        default = self.run_lame(False, EXPANDED, mu=7, lamb=40)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "expanded data"%self.domain.getDim())

    def test_lameScalar(self):
        #check default and lame assemblers agree
        default = self.run_lame(False, SCALAR)
        lame = self.run_lame(True, SCALAR)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "scalar data"%self.domain.getDim())
        #reverse order, ensure default assembler still operational
        lame = self.run_lame(True, SCALAR, mu=7, lamb=40)
        default = self.run_lame(False, SCALAR, mu=7, lamb=40)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "scalar data"%self.domain.getDim())

    def test_lameConstant(self):
        #check default and lame assemblers agree
        default = self.run_lame(False, CONSTANT)
        lame = self.run_lame(True, CONSTANT)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "constant data"%self.domain.getDim())
        #reverse order, ensure default assembler still operational
        lame = self.run_lame(True, CONSTANT, mu=7, lamb=40)
        default = self.run_lame(False, CONSTANT, mu=7, lamb=40)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "constant data"%self.domain.getDim())

class Test_RipleyCustomAssemblers2D(RipleyCustomAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(20,20)

    def tearDown(self):
        del self.domain

class Test_RipleyCustomAssemblers3D(RipleyCustomAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(10,10,10)

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

