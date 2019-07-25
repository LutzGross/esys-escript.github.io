##############################################################################
#
# Copyright (c) 2003-2019 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2019 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import *

class Test_GeneratorsOnFinley(unittest.TestCase):

    def test_test(self):
        self.assertEqual(1,1,"Oxley testing code has failed to run.")

    # This just checks that it writes to file and p4est doesn't return an error.
    # TODO: Check the file itself for errors
    def test_write_rectangle_to_vtk(self):
        domain1=Rectangle(order=2,n0=10,n1=10)
        domain1.writeToVTK("/tmp/rectangle")
        file_exists=os.path.exists("/tmp/rectangle.vtk_0000.vtu")
        self.assertEqual(file_exists,True,"Oxley rectangle.writeToVTK.")
        # Clean up
        os.remove("/tmp/rectangle.vtk_0000.vtu")
        os.remove("/tmp/rectangle.vtk.pvtu")
        os.remove("/tmp/rectangle.vtk.visit")

    # This just checks that it writes to file and p8est doesn't return an error
    # TODO: Check the file itself for errors
    def test_write_brick_to_vtk(self):
        domain1=Brick(order=2,n0=10,n1=10,n2=10)
        domain1.writeToVTK("/tmp/brick")
        file_exists=os.path.exists("/tmp/brick.vtk_0000.vtu")
        self.assertEqual(file_exists,True,"Oxley brick.writeToVTK.")
        # Clean up
        os.remove("/tmp/brick.vtk_0000.vtu")
        os.remove("/tmp/brick.vtk.pvtu")
        os.remove("/tmp/brick.vtk.visit")


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

