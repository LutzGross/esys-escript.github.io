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

# $Id: test_text.py,v 1.1 2005/01/11 04:21:53 paultcochrane Exp $

## @file test_text.py

import unittest
import sys,os,string
here = os.getcwd() + '/../../'
sys.path.append(here)
from pyvisi import *   # this should import all of the pyvisi stuff needed

from ESyS import *
import Finley

"""
Class and functions for testing the Text class
"""

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
    unittest.main()

# vim: expandtab shiftwidth=4:
