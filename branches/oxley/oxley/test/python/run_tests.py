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
import numpy as np
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import *

DONT_HAVE_BOOST_NUMPY = not hasFeature("boost_numpy")

class Test_OxleyDomain(unittest.TestCase):

    def test_test(self):
        self.assertEqual(1,1,"Oxley testing code has failed to run.")

    # This just checks that it writes to file and p4est doesn't return an error.
    def test_getDim(self):
        domain=Rectangle(order=2,n0=10,n1=10)
        self.assertEqual(domain.getDim(),2,"Oxley::Rectangle dims are wrong.")
        del domain
        domain=Brick(order=2,n0=10,n1=10,n2=10)
        self.assertEqual(domain.getDim(),3,"Oxley::Brick dims are wrong.")
        del domain

class Test_Rectangle(unittest.TestCase):

    # This just checks that it writes to file and p4est doesn't return an error.
    # TODO: Check the file for errors
    def test_write_rectangle_to_vtk(self):
        domain1=Rectangle(order=2,n0=10,n1=10)
        domain1.writeToVTK("/tmp/rectangle")
        file_exists=os.path.exists("/tmp/rectangle.vtk_0000.vtu")
        self.assertEqual(file_exists,True,"Oxley rectangle.writeToVTK.")
        # Clean up
        os.remove("/tmp/rectangle.vtk_0000.vtu")
        os.remove("/tmp/rectangle.vtk.pvtu")
        os.remove("/tmp/rectangle.vtk.visit")
        del domain1

    # def test_write_rectangle_to_vtk_tag_info(self):
    #     domain1=Rectangle(order=2,n0=10,n1=10)
    #     domain1.writeToVTK(filename="/tmp/rectangle",writeTagInfo=True)
    #     file_exists=os.path.exists("/tmp/rectangle.vtk_0000.vtu")
    #     self.assertEqual(file_exists,True,"Oxley rectangle.writeToVTK.")
    #     # Clean up
    #     os.remove("/tmp/rectangle.vtk_0000.vtu")
    #     os.remove("/tmp/rectangle.vtk.pvtu")
    #     os.remove("/tmp/rectangle.vtk.visit")
    #     del domain1

    def test_getNumVertices_Rectangle(self):
        domain1=Rectangle(order=2,n0=2,n1=2)
        num_corners = domain1.getNumVertices()
        self.assertEqual(num_corners,9,"test_getNumCorners")
        del domain1

    # @unittest.skipIf(DONT_HAVE_BOOST_NUMPY is True, "No boost numpy")
    # def test_addSurface_incorrect_input_2D(self):
    #     domain1=Rectangle(order=2,n0=3,n1=3)
    #     nx=np.array([1,2,3,4,5])
    #     ny=np.array([1,2,3,4,5,6])
    #     addSurface(domain1,nx,ny)
    #     del domain1


# class Test_Brick(unittest.TestCase):
#
#     def test_write_brick_to_vtk(self):
#         domain1=Brick(order=2,n0=10,n1=10,n2=10)
#         domain1.writeToVTK("/tmp/brick")
#         file_exists=os.path.exists("/tmp/brick.vtk_0000.vtu")
#         self.assertEqual(file_exists,True,"Oxley brick.writeToVTK.")
#         # Clean up
#         os.remove("/tmp/brick.vtk_0000.vtu")
#         os.remove("/tmp/brick.vtk.pvtu")
#         os.remove("/tmp/brick.vtk.visit")
#         del domain1
#
#     def test_write_brick_to_vtk_tag_info(self):
#         domain1=Brick(order=2,n0=10,n1=10,n2=10)
#         domain1.writeToVTK("/tmp/brick",writeTagInfo=True)
#         file_exists=os.path.exists("/tmp/brick.vtk_0000.vtu")
#         self.assertEqual(file_exists,True,"Oxley brick.writeToVTK.")
#         # Clean up
#         os.remove("/tmp/brick.vtk_0000.vtu")
#         os.remove("/tmp/brick.vtk.pvtu")
#         os.remove("/tmp/brick.vtk.visit")
#         del domain1
#
#     def test_getNumCorners_Brick(self):
#         domain1=Brick(order=2,n0=2,n1=2,n2=2)
#         num_corners = domain1.getNumVertices()
#         self.assertEqual(num_corners,27,"test_getNumCorners")
#         del domain1

    # @unittest.skipIf(DONT_HAVE_BOOST_NUMPY is True, "No boost numpy")
    # def test_addSurface_incorrect_input_2D(self):
    #     domain1=Bric Test_Brick(unittest.TestCase):
#
#     def test_write_brick_to_vtk(self):
#         domain1=Brick(order=2,n0=10,n1=10,n2=10)
#         domain1.writeToVTK("/tmp/brick")
#         file_exists=os.path.exists("/tmp/brick.vtk_0000.vtu")
#         self.assertEqual(file_exists,True,"Oxley brick.writeToVTK.")
#         # Clean up
#         os.remove("/tmp/brick.vtk_0000.vtu")
#         os.remove("/tmp/brick.vtk.pvtu")
#         os.remove("/tmp/brick.vtk.visit")
#         del domain1
#
#     def test_write_brick_to_vtk_tag_info(self):
#         domain1=Brick(order=2,n0=10,n1=10,n2=10)
#         domain1.writeToVTK("/tmp/brick",writeTagInfo=True)
#         file_exists=os.path.exists("/tmp/brick.vtk_0000.vtu")
#         self.assertEqual(file_exists,True,"Oxley brick.writeToVTK.")
#         # Clean up
#         os.remove("/tmp/brick.vtk_0000.vtu")
#         os.remove("/tmp/brick.vtk.pvtu")
#         os.remove("/tmp/brick.vtk.visit")
#         del domain1
#
#     def test_getNumCorners_Brick(self):
#         domain1=Brick(order=2,n0=2,n1=2,n2=2)
#         num_corners = domain1.getNumVertices()
#         self.assertEqual(num_corners,27,"test_getNumCorners")
#         del domain1

# def test_uniform_refine_3D(self):
#     domain1=Brick(order=2,n0=1,n1=1,n2=1)
#     initial_corners = domain1.getNumCorners();
#     domain1.refine(1,"uniform")
#     final_corners = domain1.getNumCorners()
#     self.assertEqual(initial_corners,8,"test_uniform_refine_2D")
#     self.assertEqual(initial_corners,27,"test_uniform_refine_2D")

# AEAE: List of tests to add
# addSurface works for at least one test example in 2D
# addSurface works for at least one test example in 3D
# domain.refine bounces an error when incorrect input is passed on
# uniform refinement works properly for both rectangle and bricks
# uniform refinement works properly for irregular rectangle and bricks
# writeToVTK works for a normal mesh in 2D
# writeToVTK works for a normal mesh in 3D
# writeToVTK works for a normal mesh in 2D with Tag info
# writeToVTK works for a normal mesh in 3D with Tag info


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
