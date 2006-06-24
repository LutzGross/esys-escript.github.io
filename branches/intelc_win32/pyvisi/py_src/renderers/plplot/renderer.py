"""
Brief introduction to what the file contains/does

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


from common import debugMsg, _rendererName
from esys.pyvisi.renderer import Renderer as BaseRenderer

class Renderer(BaseRenderer):
    """
    Generic object holding a renderer of a Scene()
    """

    def __init__(self):
        """
        Initialisation
        """
        debugMsg("Called Renderer.__init__()")
        BaseRenderer.__init__(self)  # initialisation of base class
    
        # initialise some attributes
        self.renderWindowWidth = 640
        self.renderWindowHeight = 480

        self.name = _rendererName

        # the namespace to run the exec code
        self.renderDict = {}

        # initialise the evalstack
        self._evalStack = ""

        # keep the initial setup of the module for later reuse
        self._initStack = ""

        # initialise the renderer module
        self.addToInitStack("# Renderer.__init__()")
        self.addToInitStack("import plplot")

        # we now need the Numeric package so can handle arrays better
        self.addToInitStack("from Numeric import *")

# vim: expandtab shiftwidth=4:
