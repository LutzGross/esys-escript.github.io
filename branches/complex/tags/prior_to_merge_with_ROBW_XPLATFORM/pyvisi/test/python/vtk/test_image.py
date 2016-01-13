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

# $Id: test_image.py,v 1.1 2005/01/11 04:21:53 paultcochrane Exp $

## @file test_image.py

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
import pyvisi   # this should import all of the pyvisi stuff needed

#from ESyS import *
#import Finley

"""
Class and functions for testing the Image classes
"""

class TestImage(unittest.TestCase):
    """
    The main test class
    """

    def testExactlyOneArg(self):
        """
        Tests that Image class is passed exactly one argument
        """
        # this tests for no args passed
        self.assertRaises(TypeError, pyvisi.Image.__init__)
        # this tests for two args passed
        self.assertRaises(TypeError, pyvisi.Image.__init__, 'one', 'two')

if __name__ == '__main__':
    unittest.main()

# vim: expandtab shiftwidth=4:
