##############################################################################
#
# Copyright (c) 2003-2019 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"
__author__="Adam Ellery, a.ellery@uq.edu.au"

import os
import numpy as np
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import *

DONT_HAVE_BOOST_NUMPY = not hasFeature("boost_numpy")

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

# class Test_LinearPDEOnOxleyRectangle(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_2Do1, Test_TransportPDE):
@unittest.skip("Oxley Rectangle meshes have SystemMatrixPattern errors with LinearPDE - see issue #118")
class Test_LinearPDEOnOxleyRectangle(Test_LinearPDE):
    RES_TOL=1.e-7
    ABS_TOL=1.e-8
    def setUp(self):
        # for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
        #     NX=x
        #     NY=mpiSize//x
        #     if NX*NY == mpiSize:
        #         break
        NX = 1
        NY = 1
        self.domain=Rectangle(n0=NE*NX-1, n1=NE*NY-1, l0=1., l1=1., d0=NX, d1=NY)
        self.order = 1
    def tearDown(self):
        del self.domain

# class Test_LinearPDEOnOxleyBrick(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_3Do1, Test_TransportPDE):
#TODO
# class Test_LinearPDEOnOxleyBrick(Test_LinearPDE, Test_LameEquation, Test_Helmholtz, Test_LinearPDE_noLumping, Test_pdetools, Test_assemblage_3Do1):
#     RES_TOL=1.e-7
#     ABS_TOL=1.e-8
#     def setUp(self):
#         for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
#             NX=x[0]
#             NY=x[1]
#             NZ=mpiSize//(x[0]*x[1])
#             if NX*NY*NZ == mpiSize:
#                 break

#         self.domain = Brick(n0=NE*NX-1, n1=NE*NY-1, n2=NE*NZ-1, l0=1., l1=1., l2=1., d0=NX, d1=NY, d2=NZ)
#         self.order = 1

#     def tearDown(self):
#         del self.domain

@unittest.skip("Oxley Rectangle meshes have SystemMatrixPattern errors with Poisson - see issue #118")
class Test_PoissonOnOxley(Test_Poisson):
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
