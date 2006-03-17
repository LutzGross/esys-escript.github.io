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

# $Id: test_renderer.py,v 1.1 2005/01/11 05:07:41 paultcochrane Exp $

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
