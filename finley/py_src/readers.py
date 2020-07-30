# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
some mesh handling

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.pycad.gmsh import Design as GMSHDesign
from .factorywrappers import ReadGmsh, ReadMesh
from .finleycpp import LoadMesh

def MakeDomain(design,integrationOrder=-1, reducedIntegrationOrder=-1, optimizeLabeling=True, useMacroElements=False):
    """
    Creates a Finley `Domain` from a `esys.pycad.design.Design` object.
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
    if isinstance(design, GMSHDesign):
        if useMacroElements: design.setElementOrder(2)
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
        raise TypeError("Finley does not support %s designs."%design.__class__.__name__)
    # fill in the tag map
    design.getTagMap().passToDomain(dom)
    return dom

def GetMeshFromFile(filename, **kwargs):
    """
    Reads a mesh from a file, determines the reader to use based on the file
    extension. All cases require a filename and gmsh files require a number
    of dimensions (it doesn't hurt to pass this in all the time). Other
    keyword args come from the underlying reader functions.
    """
    spl=filename.split('.')
    ext=spl[len(spl)-1]
    # extract possible params
    integrationOrder=-1
    if "integrationOrder" in kwargs:
        integrationOrder=kwargs["integrationOrder"]
    reducedIntegrationOrder=-1
    if "reducedIntegrationOrder" in kwargs:
        integrationOrder=kwargs["reducedIntegrationOrder"]
    optimize=True  
    if "optimize" in kwargs:
        integrationOrder=kwargs["optimize"]
    useMacroElements=False
    if "useMacroElements" in kwargs:
        integrationOrder=kwargs["useMacroElements"]
    if ext=="fly":
        return ReadMesh(filename, integrationOrder, reducedIntegrationOrder, optimize)
    elif ext=="msh":
        if not "numDim" in kwargs:
           raise ValueError("The numDim argument is required in order to read .msh files.")
        return ReadGmsh(filename, kwargs['numDim'], integrationOrder, reducedIntegrationOrder, optimize, useMacroElements)
    else:
#        return LoadMesh(filename)
        raise ValueError("Unsupported extension .%s"%ext)
