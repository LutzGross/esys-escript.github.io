
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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for the linearPDE  and pdetools test on finley

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_nonLinearPDE import Test_nlpde
from esys.finley import Rectangle, Brick

class Test_nonLinearPDEOnFinley2D(Test_nlpde):
   def setUp(self):
        self.domain = Rectangle(l0=1.,l1=1.,n0=10, n1=10) 
   def tearDown(self):
        del self.domain

class Test_nonLinearPDEOnFinley3D(Test_nlpde):
   def setUp(self):
        self.domain = Brick(l0=1.,l1=1.,l2=1.,n0=10, n1=10,n2=10) 
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

