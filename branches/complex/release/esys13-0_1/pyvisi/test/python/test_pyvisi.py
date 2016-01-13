"""
Class and functions for testing the pyvisi framework 

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


import unittest
from esys.pyvisi import *   # this should import all of the pyvisi stuff needed

class TestAxes(unittest.TestCase):
   """
   The main test class
   """

   def testAxesBases(self):
        bases = "%s" % Axes.__bases__
        self.assertEqual(bases, '<class \'esys.pyvisi.scene.Plot\'>')

class TestCamera(unittest.TestCase):
    """
    The main test class
    """

    def testAddCamera(self):
        scene = Scene()
        camera = scene.addCamera()
        self.assertEqual(camera.__class__.__name__, 'Camera')

    def testOnlyOneCamera(self):
        # there should only be one camera object
        scene = Scene()
        camera = scene.addCamera()
        self.assertEqual(camera._cameraCount,1)
        self.assertRaises(somekindoferror, scene.addCamera())


class TestImage(unittest.TestCase):
    """
    The main test class
    """

    def testExactlyOneArg(self):
        """
        Tests that Image class is passed exactly one argument
        """
        # this tests for no args passed
        self.assertRaises(TypeError, Image.__init__)
        # this tests for two args passed
        self.assertRaises(TypeError, Image.__init__, 'one', 'two')

class TestItem(unittest.TestCase):
    """
    The main test class
    """

    def testFunction(self):
        """
        A test function
        """
        self.assertEqual()
        self.assertRaises()
        self.assert_()


class TestPlane(unittest.TestCase):
    """
    The main test class
    """

    def testFunction(self):
        """
        A test function
        """
        self.assertEqual()
        self.assertRaises()
        self.assert_()



class TestPlot(unittest.TestCase):
    """
    The main test class
    """

    def testArrowPlotBases(self):
        bases = "%s" % ArrowPlot.__bases__
        self.assertEqual(bases, '<class \'esys.pyvisi.scene.Plot\'>')

    def testArrowPlotSetData(self):
        scene = Scene()
        plot = scene.addArrowPlot()
        data = Finley.Brick(3,5,7).Nodes().getX()  # old esys!!
        self.assert_(plot.setData(data),msg="Failed to set data in ArrowPlot")
    
    def testContourPlotBases(self):
        bases = "%s" % ContourPlot.__bases__
        self.assertEqual(bases, '<class \'esys.pyvisi.scene.Plot\'>')

    def testLinePlotBases(self):
        bases = "%s" % LinePlot.__bases__
        self.assertEqual(bases, '<class \'esys.pyvisi.scene.Plot\'>')


class TestRenderer(unittest.TestCase):
    """
    Test the renderer object at the base esys.pyvisi level
    """

    def setUp(self):
        self.ren = Renderer()

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


class TestScene(unittest.TestCase):
    """
    The main test class
    """

    def setUp(self):
        self.scene = Scene()

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

class TestText(unittest.TestCase):
    """
    The main test class
    """

    def testAddText(self):
        scene = Scene()
        text = scene.addText()
        self.assertEqual(text.__class__.__name__, 'Text')

    def testTextFontAttr(self):
        # make sure the font attribute exists and is text
        scene = Scene()
        text = scene.addText()
        self.assert_(text.font.isalpha(), \
            msg='Text font attribute is not text')

    def testTextChangeFont(self):
        # now try and set the font to something, and see if that works
        scene = Scene()
        text = scene.addText()
        text.font = "Helvetica"
        self.assertEqual('Helvetica', text.font)
        self.assertNotEqual('Times', text.font)

    def testTextNotAlphaNum(self):
        # text to see if when we set the font to alphanumeric that it barfs
        scene = Scene()
        text = scene.addText()
        text.font = "Times2"
        self.assert_(not text.font.isalpha(), \
            msg='Text font attribute is supposed to be text not alphanum')
    



if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(TestAxes))
   suite.addTest(unittest.makeSuite(TestCamera))
   suite.addTest(unittest.makeSuite(TestImage))
   suite.addTest(unittest.makeSuite(TestItem))
   suite.addTest(unittest.makeSuite(TestPlane))
   suite.addTest(unittest.makeSuite(TestPlot))
   suite.addTest(unittest.makeSuite(TestRenderer))
   suite.addTest(unittest.makeSuite(TestScene))
   suite.addTest(unittest.makeSuite(TestText))
   unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()
