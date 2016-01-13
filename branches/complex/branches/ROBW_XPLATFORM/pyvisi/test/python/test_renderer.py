# Copyright (C) 2004 Paul Cochrane
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# $Id: test_renderer.py,v 1.5 2005/09/20 05:57:40 paultcochrane Exp $

## @file test_plot.py

import unittest
import sys,os
from string import *
here = os.getcwd() + '/../../'
sys.path.append(here)
import pyvisi  # this should import all of the pyvisi stuff needed

"""
Class and functions for testing the Renderer class
"""

class TestRenderer(unittest.TestCase):
    """
    Test the renderer object at the base pyvisi level
    """

    def setUp(self):
        self.ren = pyvisi.Renderer()

    def testInit(self):
        """
        Tests initialisation of the Renderer object
        """
        classStr = self.ren.__class__.__name__
        self.assertEqual('Renderer', classStr)

    def testExactlyOneArgument(self):
        """
        Check that method only accepts one argument
        """
        # check setRenderWindowHeight
        self.assertRaises(TypeError, self.ren.setRenderWindowHeight)
        self.assertRaises(TypeError, self.ren.setRenderWindowHeight, 10, 20)

        # check setRenderWindowWidth
        self.assertRaises(TypeError, self.ren.setRenderWindowWidth)
        self.assertRaises(TypeError, self.ren.setRenderWindowWidth, 10, 20)

        # check addToEvalStack
        self.assertRaises(TypeError, self.ren.addToEvalStack)
        self.assertRaises(TypeError, self.ren.addToEvalStack, "moo", "baa")

        # check addToInitStack
        self.assertRaises(TypeError, self.ren.addToInitStack)
        self.assertRaises(TypeError, self.ren.addToInitStack, "moo", "baa")

    def testExactlyTwoArguments(self):
        """
        Check that method only accepts two arguments
        """
        self.assertRaises(TypeError, \
                self.ren.setRenderWindowDimensions)
        self.assertRaises(TypeError, \
                self.ren.setRenderWindowDimensions, 12)
        self.assertRaises(TypeError, \
                self.ren.setRenderWindowDimensions, 12, 14, 16)

    def testRenderWindowDefaultDims(self):
        """
        Check render window default width and height
        """
        self.assertEqual(640, self.ren.renderWindowWidth)  # width
        self.assertEqual(480, self.ren.renderWindowHeight) # height

    def testGetRenderWindowWidth(self):
        """
        Check getting the render window width
        """
        width = self.ren.getRenderWindowWidth()
        self.assertEqual(self.ren.renderWindowWidth, width)

    def testSetRenderWindowWidth(self):
        """
        Test setting the render window width
        """
        width = 720
        self.ren.setRenderWindowWidth(width)
        self.assertEqual(self.ren.getRenderWindowWidth(), width)

    def testGetRenderWindowHeight(self):
        """
        Check getting the render window height
        """
        height = self.ren.getRenderWindowHeight()
        self.assertEqual(self.ren.renderWindowHeight, height)

    def testSetRenderWindowHeight(self):
        """
        Test setting the render window height
        """
        height = 593
        self.ren.setRenderWindowHeight(height)
        self.assertEqual(self.ren.getRenderWindowHeight(), height)

    def testGetRenderWindowDimensions(self):
        """
        Test getting the render window width and height
        """
        (width, height) = self.ren.getRenderWindowDimensions()
        self.assertEqual((640,480), (width, height))

    def testSetRenderWindowDimensions(self):
        """
        Test setting the render window width and height
        """
        width = 123
        height = 456
        self.ren.setRenderWindowDimensions(width, height)
        self.assertEqual(self.ren.getRenderWindowDimensions(), (width,height))

    def testIntegerArgsToSet(self):
        """
        Test setting of integer arguments
        """
        self.assertRaises(AssertionError, \
                self.ren.setRenderWindowWidth, "moo")
        self.assertRaises(AssertionError, \
                self.ren.setRenderWindowHeight, "moo")
        self.assertRaises(AssertionError, \
                self.ren.setRenderWindowDimensions, "moo", "baa")

    def testStringArgsToSet(self):
        """
        Test setting of string arguments
        """
        self.assertRaises(AssertionError, self.ren.addToEvalStack, 10)
        self.assertRaises(AssertionError, self.ren.addToInitStack, 10)

    def testGetEvalStack(self):
        """
        Test getting the evaluation stack
        """
        # it should be the null string on initialisation
        self.assertEqual("", self.ren.getEvalStack())

    def testAddToEvalStack(self):
        """
        Test adding a string to the evaluation stack
        """
        inString = "my string"
        outString = inString + '\n'
        self.ren.addToEvalStack(inString)
        self.assertEqual(self.ren.getEvalStack(), outString)

    def testGetInitStack(self):
        """
        Test getting the initialisation stack
        """
        # should be the null string initially
        self.assertEqual("", self.ren.getInitStack())

    def testAddToInitStack(self):
        """
        Test adding a string to the initialisation stack
        """
        inString = "my string"
        outString = inString + '\n'
        self.ren.addToInitStack(inString)
        self.assertEqual(self.ren.getInitStack(), outString)

    def testResetEvalStack(self):
        """
        Test resetting the evaluation stack
        """
        # set the stack to something
        inString = "my string"
        self.ren.addToEvalStack(inString)
        # reset the stack
        self.ren.resetEvalStack()
        # now check that it's the null string again
        self.assertEqual("", self.ren.getEvalStack())

    def testResetInitStack(self):
        """
        Test resetting the initialisation stack
        """
        # set the stack to something
        inString = "my string"
        self.ren.addToInitStack(inString)
        # reset the stack
        self.ren.resetInitStack()
        # now check that it's the null string again
        self.assertEqual("", self.ren.getInitStack())

if __name__ == '__main__':
    unittest.main()

# vim: expandtab shiftwidth=4:
