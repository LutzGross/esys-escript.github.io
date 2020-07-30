
##############################################################################
#
# Copyright (c) 2014-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2014-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from .dudleycpp import __Brick_driver, __Rectangle_driver

def Rectangle(n0=1, n1=1, order=1, l0=1.0, l1=1.0, periodic0=False,
              periodic1=False, integrationOrder=-1, reducedIntegrationOrder=-1,
              useElementsOnFace=False, useFullElementOrder=False,
              optimize=False, **kwargs):
    """
    Creates a triangular mesh by subdividing n0 x n1 rectangular elements over
    the brick [0,l0] x [0,l1].
    The following keyword arguments are understood:
      diracPoints  - coordinates of dirac points to add to domain
      diracTags    - list of tags for the dirac points
      escriptworld - MPI (sub)world to use

    :param n0: number of elements for side 0
    :type n0: ``int``
    :param n1: number of elements for side 1
    :type n1: ``int``
    :param order: for compatibility with finley, always 1
    :param l0: length of side 0
    :type l0: ``float``
    :param l1: length of side 1
    :type l1: ``float``
    :param periodic0: for compatibility with finley, always False
    :param periodic1: for compatibility with finley, always False
    :param integrationOrder: for compatibility with finley, always 2
    :param reducedIntegrationOrder: for compatibility with finley, unused
    :param useElementsOnFace:  for compatiblity with finley, always False
    :param useFullElementOrder: for compatibility with finley, always False
    :param optimize: Enable optimisation of node labels
    :type optimize: ``bool``
    """
    if 'diracPoints' in kwargs:
        points=kwargs['diracPoints']
    if 'diracTags' in kwargs:
        tags=kwargs['diracTags']
    faceon=useElementsOnFace
    args=[n0, n1, order, l0, l1, periodic0, periodic1, integrationOrder, 
      reducedIntegrationOrder, faceon, useFullElementOrder, optimize];
    if 'escriptworld' in kwargs:
      args+=[kwargs['escriptworld']]
    else:
      args+=[None]
    return __Rectangle_driver(args)

Rectangle.__doc__=__Rectangle_driver.__doc__

def Brick(n0=1, n1=1, n2=1, order=1, l0=1.0, l1=1.0, l2=1.0, periodic0=False,
          periodic1=False, periodic2=False, integrationOrder=-1,
          reducedIntegrationOrder=-1, useElementsOnFace=False,
          useFullElementOrder=False, optimize=False, **kwargs):
    """
    Creates a tetrahedral mesh by subdividing n0 x n1 x n2 rectangular elements
    over the brick [0,l0] x [0,l1] x [0,l2].
    The following keyword arguments are understood:
      diracPoints  - coordinates of dirac points to add to domain
      diracTags    - list of tags for the dirac points
      escriptworld - MPI (sub)world to use

    :param n0: number of elements for side 0
    :type n0: ``int``
    :param n1: number of elements for side 1
    :type n1: ``int``
    :param n2: number of elements for side 2
    :type n2: ``int``
    :param order: for compatibility with finley, always 1
    :param l0: length of side 0
    :type l0: ``float``
    :param l1: length of side 1
    :type l1: ``float``
    :param l2: length of side 2
    :type l2: ``float``
    :param periodic0: for compatibility with finley, always False
    :param periodic1: for compatibility with finley, always False
    :param periodic2: for compatibility with finley, always False
    :param integrationOrder: for compatibility with finley, always 2
    :param reducedIntegrationOrder: for compatibility with finley, unused
    :param useElementsOnFace:  for compatiblity with finley, always False
    :param useFullElementOrder: for compatibility with finley, always False
    :param optimize: Enable optimisation of node labels
    :type optimize: ``bool``
    """
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

Brick.__doc__=__Brick_driver.__doc__

