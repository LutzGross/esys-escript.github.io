"""
Class and functions for testing the Renderer class

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
import sys,os
from string import *
here = os.getcwd() + '/../../'
sys.path.append(here)
import pyvisi  # this should import all of the pyvisi stuff needed


class TestRenderer(unittest.TestCase):
    """
    The main test class
    """

    def testRendererExactlyTwoArgs(self):
        """
        Tests that the Renderer object is instantiated only with two arguments
        """
        # test just the one argument
        self.assertRaises(TypeError, \
                pyvisi.Renderer.__init__, 'one')
        # test three arguments
        self.assertRaises(TypeError, \
                pyvisi.Renderer.__init__, 'one', 'two', 'three')

    def testRendererNameType(self):
        """
        Tests the type of the renderer name; should be alphanumeric
        """
        ren = pyvisi.Renderer()
        renName = ren.rendererName
        self.assert_(renName.isalnum(),\
                msg='Renderer() argument is not alphanumeric')

    def testRendererReturn(self):
        """
        Tests that a Renderer() object is returned
        """
        ren = pyvisi.Renderer()
        classStr = ren.__class__.__name__
        self.assertEqual('Renderer', classStr)

if __name__ == '__main__':
    unittest.main()

# vim: expandtab shiftwidth=4:
