
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

import os, sys
import esys.escriptcore.utestselect as unittest
from esys.escript import *
from esys.ripley import Rectangle, Brick
from esys.escript.linearPDEs import LameEquation

try:
     RIPLEY_WORKDIR=os.environ['RIPLEY_WORKDIR']
except KeyError:
     RIPLEY_WORKDIR='.'

EXPANDED, SCALAR, CONSTANT = range(3)

def run_lame(domain, fast, test_type, mu=3, lamb=50):
    d=domain.getDim()
    mypde=LameEquation(domain, useFast=fast)
    cf=ContinuousFunction(domain)
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
    mypde.setValue(q=msk,r=u_ex,lame_mu=mu,lame_lambda=lamb,f=(2*3+50*d)*FunctionOnBoundary(domain).getNormal())
    return mypde.getSolution()

class Test_CustomAssemblers(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_lame(self):
        labels = ["expanded", "scalar", "constant"]
        for test_type in [EXPANDED, SCALAR, CONSTANT]:
            for dom, dim in [(Brick(10,10,10), "3D"), (Rectangle(20,20), "2D")]:
                #check default and lame assemblers agree
                dom = Brick(10,10,10)
                default = run_lame(dom, False, test_type)
                lame = run_lame(dom, True, test_type)
                self.assertLess(Lsup(default - lame), 1e-8,
                        "Default and Lame " + dim + "assembler solutions " \
                        "differ for %s data"%labels[test_type])
                #reverse order, ensure default assembler still operational
                lame = run_lame(dom, True, test_type, mu=7, lamb=40)
                default = run_lame(dom, False, test_type, mu=7, lamb=40)
                self.assertLess(Lsup(default - lame), 1e-8,
                        "Default and Lame " + dim + "assembler solutions " \
                        "differ for %s data"%labels[test_type])


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_CustomAssemblers))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful():
        sys.exit(1)

