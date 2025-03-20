# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
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

from .factorywrappers import ReadGmsh, ReadMesh
from .finleycpp import LoadMesh


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
    points = []
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    tags = [ ]
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']

    if ext=="fly":
        return ReadMesh(filename, integrationOrder, reducedIntegrationOrder, optimize,  diracPoints=points, diracTags=tags)
    elif ext=="msh":
        if not "numDim" in kwargs:
           raise ValueError("The numDim argument is required in order to read .msh files.")
        return ReadGmsh(filename, kwargs['numDim'], integrationOrder, reducedIntegrationOrder, optimize, useMacroElements, diracPoints=points, diracTags=tags)
    else:
#        return LoadMesh(filename)
        raise ValueError("Unsupported extension .%s"%ext)
