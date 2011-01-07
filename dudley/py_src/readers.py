
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
some mesh handling

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au, Joel Fenwick"

from esys.pycad.gmsh import Design as GMSHDesign
from dudleycpp import ReadGmsh

def MakeDomain(design,integrationOrder=-1, reducedIntegrationOrder=-1, optimizeLabeling=True, useMacroElements=False):
    """
    Creates a Dudley `Domain` from a `esys.pycad.design.Design` object.
    Currently only gmsh is supported.

    :param design: the geometry
    :type design: `esys.pycad.design.Design`
    :param integrationOrder: integration order. If -1 the default is used.
    :type integrationOrder: ``int``
    :param reducedIntegrationOrder: reduced integration order. If -1 the
                                    default is used.
    :type reducedIntegrationOrder: ``int``
    :param optimizeLabeling: if set the labeling of the mesh nodes is optimized
    :type optimizeLabeling: ``bool``
    :param useMacroElements: uses macro elements.
    :type useMacroElements: ``bool``
    :return: the Finley domain defined by the design
    :rtype: `Domain`
    """
    if useMacroElements:
	raise TypeError("Dudley does not support macro elements")
    if isinstance(design, GMSHDesign):
        design.setElementOrder(1)
        ff=design.getFileFormat()
        design.setFileFormat(design.GMSH)
        mshname=design.getMeshHandler()
        dom = ReadGmsh(mshname,
                       design.getDim(),
                       integrationOrder,
                       reducedIntegrationOrder,
                       optimizeLabeling,
                       useMacroElements)
        design.setFileFormat(ff)
    else:
        raise TypeError("Dudley does not support %s designs."%design.__class__.__name__)
    # fill in the tag map
    design.getTagMap().passToDomain(dom)
    return dom

