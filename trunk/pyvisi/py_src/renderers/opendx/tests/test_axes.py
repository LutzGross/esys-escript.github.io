# Copyright (C) 2004-2005 Paul Cochrane
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

# $Id: test_axes.py,v 1.1 2005/02/08 05:57:56 paultcochrane Exp $

## @file test_axes.py

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
# this should import all of the pyvisi stuff needed
from pyvisi import *   
# this should import the renderer specific stuff
from pyvisi.renderers.opendx import * 
    
"""
Class and functions for testing the Axes class
"""

class TestAxes(unittest.TestCase):
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

if __name__ == '__main__':
    unittest.main()

    
# vim: expandtab shiftwidth=4:
