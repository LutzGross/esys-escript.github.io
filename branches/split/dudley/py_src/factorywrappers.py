
##############################################################################
#
# Copyright (c) 2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# This file copied and adapted from the equivalent factorywrappers.py in finley

from .dudleycpp import __Brick_driver, __Rectangle_driver

def Rectangle(n0=1, n1=1, order=1, l0=1.0, l1=1.0, periodic0=False, periodic1=False, integrationOrder=-1, 
      reducedIntegrationOrder=-1, useElementsOnFace=0, useFullElementOrder=0, optimize=0, **kwargs):
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    faceon=useElementsOnFace
    args=[n0, n1, order, l0, l1, periodic0, periodic1, integrationOrder, 
      reducedIntegrationOrder, faceon, useFullElementOrder, optimize];
    if 'escriptworld' in kwargs:
      print kwargs
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __Rectangle_driver(args)

def Brick(n0=1, n1=1, n2=1, order=1, l0=1.0, l1=1.0, l2=1.0, periodic0=0, periodic1=0, periodic2=0,
    integrationOrder=-1, reducedIntegrationOrder=-1, useElementsOnFace=0, useFullElementOrder=0,
    optimize=0, **kwargs):
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    faceon=useElementsOnFace
    args=[n0, n1, n2, order, l0, l1, l2, periodic0,  periodic1, periodic2,
    integrationOrder, reducedIntegrationOrder, faceon, useFullElementOrder,
    optimize];
    if 'escriptworld' in kwargs:
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __Brick_driver(args)
