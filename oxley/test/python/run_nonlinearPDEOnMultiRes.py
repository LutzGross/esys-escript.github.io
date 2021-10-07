
########################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for nonlinearPDEs on Oxley
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_nonLinearPDE import Test_nlpde
from esys.escript import getMPISizeWorld
# from esys.oxley import MultiResolutionDomain
from esys.oxley import Rectangle, Brick

mpiSize = getMPISizeWorld()

def test_Rectangle(**kwargs):
    m = Rectangle(**kwargs)
    m.setRefinementLevel(1)
    m.refineMesh("uniform")
    return m

# def Brick(**kwargs):
#     m = MultiResolutionDomain(3, **kwargs)
#     return m.getLevel(1)

class Test_OxleyNonLinearPDE2D(Test_nlpde):
   def setUp(self):
        self.domain = test_Rectangle(l0=1.,l1=1., n0=10, n1=10*getMPISizeWorld()-1, d1=getMPISizeWorld()) 
   def tearDown(self):
        del self.domain

# TODO
# @unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
# class Test_OxleyNonLinearPDE3D(Test_nlpde):
#    def setUp(self):
#         self.domain = Brick(l0=1.,l1=1.,l2=1., n0=10, n1=10*getMPISizeWorld()-1, n2=10, d1=getMPISizeWorld()) 
#    def tearDown(self):
#         del self.domain

if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)

