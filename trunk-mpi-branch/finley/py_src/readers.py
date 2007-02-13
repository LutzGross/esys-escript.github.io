# $Id$
"""
some mesh handling 

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006, 2007 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from esys.escript import *
from esys.pycad.gmsh import Design as GMSHDesign
from finleycpp import ReadGmsh

def MakeDomain(design,integrationOrder=-1, reducedIntegrationOrder=-1, optimizeLabeling=True):
        """
        creates a Finley L{Domain} from a L{esys.pycad.design.Design} object. Currently only
        gmsh is supported.

        @param design: the geometry
        @type design: L{esys.pycad.design.Design}
        @param integrationOrder: integration order. If -1 the default is used.
        @type integrationOrder: C{int}
        @param reducedIntegrationOrder: reduced integration order. If -1 the default is used.
        @type reducedIntegrationOrder: C{int}
        @param optimizeLabeling: if set the labeling of the mesh nodes is optimized
        @type: C{bool}
        @return: the Finley domain defined by the designs 
        @rtype: L{Domain}
        """
        if isinstance(design, GMSHDesign):
            mshname=design.getMeshHandler()
            dom = ReadGmsh(mshname,
                           design.getDim(),
                           integrationOrder,
                           reducedIntegrationOrder,
                           optimizeLabeling)
        else: 
            raise TypeError("Finley does not support %s designs."%design.__class__.__name__)
        return dom
