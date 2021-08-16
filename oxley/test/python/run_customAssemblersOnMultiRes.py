
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle, Brick
from esys.escript.linearPDEs import LameEquation

from run_customAssemblersOnOxley import OxleyLameAssemblerTestBase, OxleyWaveAssemblerTestBase, Ricker

mpiSize = getMPISizeWorld()

# def Rectangle(**kwargs):
#     m = Rectangle(n0=10,n1=10,l0=10,l1=10)
#     m.refineRegion(x0=3,x1=7,y0=4,y1=8)
#     return m

# def Brick(**kwargs):
#     m = Brick(n0=10,n1=10,n2=10,l0=10,l1=10,l2=10)
#     m.refineRegion(x0=3,x1=7,x2=7,y0=4,y1=8,y2=6)
#     return m

class Test_OxleyWaveAssembler2D(OxleyWaveAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(n0=20,n1=20,l0=100.,l1=100., diracTags=["source"], diracPoints=[(0,0)])
        self.wavelet = Ricker(100.)

        
    def tearDown(self):
        del self.domain

@unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
class Test_OxleyWaveAssembler3D(OxleyWaveAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(n0=10,n1=10,n2=10,l0=100.,l1=100., l2=100., diracTags=["source"], diracPoints=[(0,0,0)])
        self.wavelet = Ricker(100.)

    def tearDown(self):
        del self.domain

class Test_OxleyLameAssemblers2D(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(n0=20,n1=20)

    def tearDown(self):
        del self.domain

@unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
class Test_OxleyLameAssemblers3D(OxleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(n0=10,n1=10,n2=10)

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

