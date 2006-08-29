# $Id$
"""
Test suite for the pyvisi

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

import unittest
from pyvisi import *

class Test_pyvisi(unittest.TestCase):
    def test_camera(self):
      c=Camera(self.scene)
      self.failUnless(abs(c.getPosition()[0])=0.,"wrong default camera position")

