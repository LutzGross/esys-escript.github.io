
########################################################
#
# Copyright (c) 2003-2023 by The University of Queensland
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
__url__="https://launchpad.net/escript-finley"

"""
Test suite for the RefinementZone class

:remark:

"""

__author__="Adam Ellery, a.ellery@uq.edu.au"

# import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle, Brick, RefinementZone2D, RefinementZone3D

N0=10
N1=10
N2=10
L0=10.
L1=10.
L2=10.

class Test_TableRefinementZone(unittest.TestCase):
    @unittest.skip("Oxley RefinementZone causes null pointer in JMPI shared_ptr - see issue #118")
    def test_RefinementZone_2D(self):
    	domain=Rectangle(n0=N0,n1=N1,l0=L0,l1=L1)
    	zone=RefinementZone2D()

    	domain.setRefinementLevel(4)
    	domain.refinePoint(x0=0.5,y0=0.5)
    	zone.refinePoint(x0=0.5,y0=0.5,level=4)

    	domain.setRefinementLevel(3)
    	domain.refineRegion(x0=3,y0=3)
    	zone.refinePoint(x0=3,y0=3,level=3)

    	domain2=domain.applyRefinementZone(zone)

    	self.assertTrue(domain==domain2)

    # def test_RefinementZone_3D(self):
    # 	domain=Brick(n0=N0,n1=N1,n2=N2,l0=L0,l1=L1,l2=L2)
    # 	zone=RefinementZone2D()

    # 	domain.setRefinementLevel(4)
    # 	domain.refinePoint(x0=0.5,y0=0.5,z0=0.5)
    # 	zone.refinePoint(x0=0.5,y0=0.5,z0=0.5,level=4)

    # 	domain.setRefinementLevel(3)
    # 	domain.refineRegion(x0=3,y0=3,z0=3,x1=6,y1=6,z1=6)
    # 	zone.refinePoint(x0=3,y0=3,z0=3,x1=6,y1=6,z1=6,level=3)

    # 	domain2=domain.applyRefinementZone(zone)

    # 	self.assertTrue(domain==domain2)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)