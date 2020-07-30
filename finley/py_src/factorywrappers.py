
##############################################################################
#
# Copyright (c) 2011-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2011-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


from .finleycpp import __Brick_driver, __Rectangle_driver, __ReadMesh_driver, __ReadGmsh_driver


def ReadMesh(filename, integrationOrder=-1, reducedIntegrationOrder=-1, optimize=True, **kwargs):
    points=[]
    tags=[]
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    args=[filename, integrationOrder, reducedIntegrationOrder, optimize, points, tags];
    if 'escriptworld' in kwargs:
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __ReadMesh_driver(args)
  
ReadMesh.__doc__=__ReadMesh_driver.__doc__  
  
def ReadGmsh(fileName, numDim, integrationOrder=-1, reducedIntegrationOrder=-1, optimize=True,  
      useMacroElements=False, **kwargs):
    points=[]
    tags=[]
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    args=[fileName, numDim, integrationOrder, reducedIntegrationOrder, optimize,  
      useMacroElements, points, tags];
    if 'escriptworld' in kwargs:
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __ReadGmsh_driver(args)      

ReadGmsh.__doc__=__ReadGmsh_driver.__doc__


def Rectangle(n0=1, n1=1, order=1, l0=1.0, l1=1.0, periodic0=False, periodic1=False, integrationOrder=-1, 
      reducedIntegrationOrder=-1, useElementsOnFace=None, useFullElementOrder=0, optimize=0, **kwargs):
    points=[]
    tags=[]
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    faceon=useElementsOnFace
    if useElementsOnFace is None:       #We want to use 1 as the default, but only where it makes sense
        if useFullElementOrder or order==-1:
            faceon=0    #Don't use it
        else:
            faceon=1
    args=[n0, n1, order, l0, l1, periodic0, periodic1, integrationOrder, 
      reducedIntegrationOrder, faceon, useFullElementOrder, optimize, points, tags];
    if 'escriptworld' in kwargs:
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __Rectangle_driver(args)

Rectangle.__doc__=__Rectangle_driver.__doc__

def Brick(n0=1, n1=1, n2=1, order=1, l0=1.0, l1=1.0, l2=1.0, periodic0=0, periodic1=0, periodic2=0,
    integrationOrder=-1, reducedIntegrationOrder=-1, useElementsOnFace=1, useFullElementOrder=0,
    optimize=0, **kwargs):
    points=[]
    tags=[]
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    faceon=useElementsOnFace
    if useElementsOnFace is None:       #We want to use 1 as the default, but only where it makes sense
        if useFullElementOrder or order==-1:
            faceon=0    #Don't use it
        else:
            faceon=1
    args=[n0, n1, n2, order, l0, l1, l2, periodic0,  periodic1, periodic2,
    integrationOrder, reducedIntegrationOrder, faceon, useFullElementOrder,
    optimize, points, tags];
    if 'escriptworld' in kwargs:
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __Brick_driver(args)

Brick.__doc__=__Brick_driver.__doc__
