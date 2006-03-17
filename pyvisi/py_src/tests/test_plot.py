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

# $Id: test_plot.py,v 1.2 2004/11/24 06:24:06 paultcochrane Exp $

## @file test_plot.py

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
from pyvisi import *   # this should import all of the pyvisi stuff needed

from ESyS import *
import Finley

"""
Class and functions for testing the Plot classes
"""

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
