"""
Class and functions for testing the Plot classes

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

import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
from pyvisi import *   # this should import all of the pyvisi stuff needed

from ESyS import *
import Finley


class TestPlot(unittest.TestCase):
    """
    The main test class
    """

    def testArrowPlotBases(self):
        bases = "%s" % ArrowPlot.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

    def testArrowPlotSetData(self):
        scene = Scene()
        plot = scene.addArrowPlot()
        data = Finley.Brick(3,5,7).Nodes().getX()  # old esys!!
        self.assert_(plot.setData(data),msg="Failed to set data in ArrowPlot")
    
    def testContourPlotBases(self):
        bases = "%s" % ContourPlot.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

    def testLinePlotBases(self):
        bases = "%s" % LinePlot.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

if __name__ == '__main__':
    unittest.main()

# vim: expandtab shiftwidth=4:
