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

# $Id: test_scene.py,v 1.3 2005/09/20 06:41:44 paultcochrane Exp $

## @file test_scene.py

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
import pyvisi   # this should import all of the pyvisi stuff needed

"""
Class and functions for testing the Scene class
"""

class TestScene(unittest.TestCase):
    """
    The main test class
    """

    def setUp(self):
        self.scene = pyvisi.Scene()

    def testInit(self):
        """
        Check initialisation
        """
        # check that we get a Scene object
        classStr = self.scene.__class__.__name__
        self.assertEqual('Scene', classStr)

        # check default values of xSize and ySize
        self.assertEqual(640, self.scene.xSize)
        self.assertEqual(480, self.scene.ySize)

        # check that we get a renderer
        classStr = self.scene.renderer.__class__.__name__
        self.assertEqual('Renderer', classStr)

    def testAdd(self):
        """
        Test adding an item to the scene
        """
        # if no object passed, should barf with TypeError, or ValueError
        self.assertRaises((TypeError,ValueError), self.scene.add)

        # calling this at the base level should fail
        self.assertRaises(NotImplementedError, self.scene.add, "")

    def testDelete(self):
        """
        Test deleting an item from the scene
        """
        # if no object passed, should barf with TypeError, or ValueError
        self.assertRaises((TypeError,ValueError), self.scene.delete)

        # calling this at the base level should fail
        self.assertRaises(NotImplementedError, self.scene.delete, "")

    def testPlace(self):
        """
        Test placing an item within the scene
        """
        # if no object passed, should barf with TypeError, or ValueError
        self.assertRaises((TypeError,ValueError), self.scene.place)

        # calling this at the base level should fail
        self.assertRaises(NotImplementedError, self.scene.place, "")

    def testRenderArguments(self):
        """
        Test the arguments passed to the render method
        """
        self.assertRaises(ValueError, self.scene.render, pause="")
        self.assertRaises(ValueError, self.scene.render, interactive="")

    def testSaveArguments(self):
        """
        Test the arguments passed to the save method
        """
        self.assertRaises(TypeError, self.scene.save)
        self.assertRaises(ValueError, self.scene.save, fname="", format="")
        self.assertRaises((TypeError, ValueError), self.scene.save, fname="")
        self.assertRaises((TypeError, ValueError), self.scene.save, format="")
        self.assertRaises((TypeError, ValueError), \
                self.scene.save, fname="moo.file")
        self.assertRaises((TypeError, ValueError), \
                self.scene.save, format="myformat")

        # need to barf if this method is used from the base pyvisi class
        self.assertRaises(NotImplementedError, \
                self.scene.save, fname="moo.file", format="myformat")

    def testSetBackgroundColor(self):
        """
        Test setting the scene's background colour
        """
        self.assertRaises(NotImplementedError, self.scene.setBackgroundColor)

    def testGetBackgroundColor(self):
        """
        Test getting the scene's background colour
        """
        self.assertRaises(NotImplementedError, self.scene.getBackgroundColor)

    def testGetSize(self):
        """
        Test getting the size of the scene dimensions
        """
        self.assertEqual((640,480), self.scene.getSize())

    def testSetSize(self):
        """
        Test setting the size of the scene dimensions
        """
        self.scene.setSize(123, 456)
        self.assertEqual((123, 456), self.scene.getSize())

    def testRendererCommand(self):
        """
        Test the renderer command
        """
        inCommand = "mooo"
        outCommand = inCommand + '\n'
        self.scene.rendererCommand(inCommand)
        self.assertEqual(self.scene.renderer.getEvalStack(), outCommand)

    def testExactlyNoArguments(self):
        """
        Test methods that do not take an argument at all
        """
        self.assertRaises(TypeError, self.scene.getBackgroundColor, "")
        self.assertRaises(TypeError, self.scene.getSize, "")

    def testExactlyOneArgument(self):
        """
        Test methods that take exactly one argument
        """
        self.assertRaises(TypeError, self.scene.add)
        self.assertRaises(TypeError, self.scene.add, "", "")
        self.assertRaises(TypeError, self.scene.delete)
        self.assertRaises(TypeError, self.scene.delete, "", "")
        self.assertRaises(TypeError, self.scene.place)
        self.assertRaises(TypeError, self.scene.place, "", "")
        self.assertRaises(TypeError, self.scene.rendererCommand)
        self.assertRaises(TypeError, self.scene.rendererCommand, "", "")

    def testExactlyTwoArguments(self):
        """
        Test methods that take exactly two arguments
        """
        # scene.setSize()
        self.assertRaises(TypeError, self.scene.setSize)
        self.assertRaises(TypeError, self.scene.setSize, 1)
        self.assertRaises(TypeError, self.scene.setSize, 1, 2, 3)

        # scene.save()
        self.assertRaises(TypeError, self.scene.save)
        self.assertRaises(TypeError, self.scene.save, "1")
        self.assertRaises(TypeError, self.scene.save, "1", "2", "3")

    def testIntegerArguments(self):
        """
        Test methods that take integer arguments
        """
        self.assertRaises(AssertionError, self.scene.setSize, "moo", "baa")

    def testStringArguments(self):
        """
        Test methods that take string arguments
        """
        self.assertRaises(AssertionError, self.scene.save, 1, 2)
        self.assertRaises(AssertionError, self.scene.rendererCommand, 1)

if __name__ == '__main__':
    unittest.main()

# vim: expandtab shiftwidth=4:
