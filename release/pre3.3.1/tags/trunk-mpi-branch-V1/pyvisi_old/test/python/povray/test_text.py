"""
Class and functions for testing the Text class

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
