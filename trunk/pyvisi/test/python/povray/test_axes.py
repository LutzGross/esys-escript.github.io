"""
Class and functions for testing the Axes class

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

class TestAxes(unittest.TestCase):
    """
    The main test class
    """

   def testAxesBases(self):
        bases = "%s" % Axes.__bases__
        self.assertEqual(bases, '<class \'pyvisi.scene.Plot\'>')

if __name__ == '__main__':
    unittest.main()

# vim: expandtab shiftwidth=4:
