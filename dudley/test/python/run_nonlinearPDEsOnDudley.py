
########################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_nonLinearPDE import Test_nlpde
from esys.dudley import Rectangle, Brick


class Test_nonLinearPDEOnDudley2D(Test_nlpde):
   def setUp(self):
        self.domain = Rectangle(l0=1.,l1=1.,n0=10, n1=10) 
   def tearDown(self):
        del self.domain

class Test_nonLinearPDEonDudley3D(Test_nlpde):
   def setUp(self):
        self.domain = Brick(l0=1.,l1=1.,l2=1.,n0=10, n1=10, n2=10) 
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

