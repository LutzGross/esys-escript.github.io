##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
Earth Systems Science Computational Center (ESSCC)
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"
__author__="Adam Ellery, a.ellery@uq.edu.au"

"""
Test suite for the linearPDE and pdetools on oxley

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
# from test_linearPDEs import Test_Poisson, Test_LinearPDE, Test_TransportPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping
from test_linearPDEs import Test_Poisson, Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
from esys.oxley import Rectangle, Brick


try:
     OXLEY_TEST_DATA=os.environ['OXLEY_TEST_DATA']
except KeyError:
     OXLEY_TEST_DATA='.'

NE=8 # initial number of elements in each spatial direction (must be even)
mpiSize=getMPISizeWorld()

# Domain decomposition across ranks is handled internally by p4est, so the
# meshes are built from a block grid (n0/n1[/n2] blocks) regardless of mpiSize.
class Test_LinearPDEOnOxleyRectangle(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        self.domain=Rectangle(n0=NE, n1=NE, l0=1., l1=1.)
        self.order = 1
    def tearDown(self):
        del self.domain

class Test_LinearPDEOnOxleyBrick(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_3Do1):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        self.domain = Brick(n0=NE, n1=NE, n2=NE, l0=1., l1=1., l2=1.)
        self.order = 1
    def tearDown(self):
        del self.domain

class Test_PoissonOnOxleyRectangle(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        self.domain=Rectangle(n0=NE, n1=NE, l0=1., l1=1.)
    def tearDown(self):
        del self.domain

class Test_PoissonOnOxleyBrick(Test_Poisson):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        self.domain=Brick(n0=NE, n1=NE, n2=NE, l0=1., l1=1., l2=1.)
    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
