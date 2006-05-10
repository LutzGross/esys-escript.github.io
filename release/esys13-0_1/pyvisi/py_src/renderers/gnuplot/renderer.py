"""
This is the file for the Renderer class

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


# generic imports
from common import debugMsg
from common import _rendererName
from esys.pyvisi.renderer import Renderer as BaseRenderer

class Renderer(BaseRenderer):
    """
    A generic object holding a renderer of a Scene().
    """

    def __init__(self):
        """
        Initialisation of Renderer() class
        """
        debugMsg("Called Renderer.__init__()")
        BaseRenderer.__init__(self)

        # initialise some attributes
        self.renderWindowWidth = 640
        self.renderWindowHeight = 480

        # what is the name of my renderer?
        self.name = _rendererName

        # initialise the evalstack
        self._evalStack = ""

        # the namespace to run the exec code
        self.renderDict = {}

        # keep the initial setup of the module for later reuse
        self._initStack = ""

        # initialise the renderer module
        self.addToInitStack("# Renderer.__init__()")
        self.addToInitStack("import Gnuplot")
        # Gnuplot package needs Numeric package
        self.addToInitStack("from Numeric import *")

        # need to add a check here to see if the Numeric module has been
        # imported, and if not, throw an error. 

        self.addToInitStack("_gnuplot = Gnuplot.Gnuplot()")

# vim: expandtab shiftwidth=4:
