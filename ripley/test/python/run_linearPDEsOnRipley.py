
########################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# Earth Systems Science Computational Center (ESSCC)
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
########################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
Earth Systems Science Computational Center (ESSCC)
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

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
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_linearPDEs import Test_Poisson, Test_LinearPDE, Test_TransportPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_3Do1
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
from esys.ripley import Rectangle, Brick


try:
     RIPLEY_TEST_DATA=os.environ['RIPLEY_TEST_DATA']
except KeyError:
     RIPLEY_TEST_DATA='.'

NE=8 # number of element in each spatial direction (must be even)
mpiSize=getMPISizeWorld()

class Test_LinearPDEOnRipleyRect(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1, Test_TransportPDE):
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

class Test_LinearPDEOnRipleyBrick(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_3Do1, Test_TransportPDE):
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
    run_tests(__name__, exit_on_failure=True)

