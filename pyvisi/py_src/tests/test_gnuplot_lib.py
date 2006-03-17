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

# $Id: test_gnuplot_lib.py,v 1.2 2004/11/24 06:24:06 paultcochrane Exp $

## @file test_gnuplot_lib.py

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
from pyvisi import *    # this should import the general stuff from pyvisi
# shouldn't that line be from pyvisi import *??
from pyvisi.lib.gnuplot import *  # this should import the gnuplot interface

from ESyS import *
import Finley

"""
Class and functions for testing the Gnuplot interface
"""

class TestGnuplot(unittest.TestCase):
    """
    The main test class
    """

    # tests of plots
    def testAddArrowPlot(self):
        scene = Scene()
        plot = scene.addArrowPlot()
        self.assertEqual(plot.__class__.__name__, 'ArrowPlot')

    def testArrowPlotBases(self):
        bases = "%s" % ArrowPlot.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

    def testArrowPlotSetData(self):
        scene = Scene()
        plot = scene.addArrowPlot()
        data = Finley.Brick(3,5,7).Nodes().getX()  # old esys!!
        self.assert_(plot.setData(data),msg="Failed to set data in ArrowPlot")
    
    def testAddContourPlot(self):
        scene = Scene()
        plot = scene.addContourPlot()
        self.assertEqual(plot.__class__.__name__, 'ContourPlot')

    def testContourPlotBases(self):
        bases = "%s" % ContourPlot.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

    def testAddLinePlot(self):
        scene = Scene()
        plot = scene.addLinePlot()
        self.assertEqual(plot.__class__.__name__, 'LinePlot')

    def testLinePlotBases(self):
        bases = "%s" % LinePlot.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

    # end of tests of plots

    def testAddAxes(self):
        scene = Scene()
        axes = scene.addAxes()
        self.assertEqual(axes.__class__.__name__, 'Axes')

    def testAxesBases(self):
        bases = "%s" % Axes.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

    def testAddText(self):
        scene = Scene()
        text = scene.addText()
        self.assertEqual(text.__class__.__name__, 'Text')

    def testAddCamera(self):
        scene = Scene()
        camera = scene.addCamera()
        self.assertEqual(camera.__class__.__name__, 'Camera')

    def testAddImage(self):
        scene = Scene()
        image = scene.addImage()
        self.assertEqual(image.__class__.__name__, 'Image')

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
    unittest.main()

# vim: expandtab shiftwidth=4:
