# $Id$
"""
runs the pyvisi tests

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

from test_pyvisi import *

class Test_vtk_jpeg(Test_pyvisi):
    def setUp(self):
        self.scene=OpenScene(renderer="vtk_jpeg")

class Test_vtk_ps(Test_pyvisi):
    def setUp(self):
        self.scene=OpenScene(renderer="vtk_ps")



if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_vtk_jpeg))
   suite.addTest(unittest.makeSuite(Test_vtk_ps))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if s.wasSuccessful():
     sys.exit(0)
   else:
     sys.exit(1)
