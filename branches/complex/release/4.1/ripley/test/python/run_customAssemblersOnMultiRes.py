
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.ripley import MultiResolutionDomain
from esys.escript.linearPDEs import LameEquation

from run_customAssemblersOnRipley import RipleyLameAssemblerTestBase, RipleyWaveAssemblerTestBase, Ricker

mpiSize = getMPISizeWorld()

def Rectangle(**kwargs):
    m = MultiResolutionDomain(2, **kwargs)
    return m.getLevel(1)

def Brick(**kwargs):
    m = MultiResolutionDomain(3, **kwargs)
    return m.getLevel(1)

class Test_RipleyWaveAssembler2D(RipleyWaveAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(n0=20,n1=20,l0=100.,l1=100., diracTags=["source"],
                diracPoints=[(0,0)])
        self.wavelet = Ricker(100.)

        
    def tearDown(self):
        del self.domain

@unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
class Test_RipleyWaveAssembler3D(RipleyWaveAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(n0=10,n1=10,n2=10,l0=100.,l1=100., l2=100.,
                diracTags=["source"], diracPoints=[(0,0,0)])
        self.wavelet = Ricker(100.)

    def tearDown(self):
        del self.domain


class Test_RipleyLameAssemblers2D(RipleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(n0=20,n1=20)

    def tearDown(self):
        del self.domain

@unittest.skipIf(mpiSize > 1, "3D Multiresolution domains require single process")
class Test_RipleyLameAssemblers3D(RipleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(n0=10,n1=10,n2=10)

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

