
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
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################


__copyright__="""Copyright (c) 2011-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


from .finleycpp import __Brick_driver, __Brick_driver_MPI, __Rectangle_driver, __Rectangle_driver_MPI, __ReadMesh_driver, __ReadGmsh_driver


def ReadMesh(filename, integrationOrder=-1, reducedIntegrationOrder=-1, optimize=True,
             diracPoints=[], diracTags=[], comm=None):
    """
    Read a mesh from a file. For MPI parallel runs fan out the mesh to multiple processes.

    :param filename: Path to mesh file
    :param integrationOrder: Order of quadrature scheme (-1 for automatic)
    :param reducedIntegrationOrder: Order of reduced quadrature scheme (-1 for automatic)
    :param optimize: Enable optimization of node labels
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, defaults to MPI_COMM_WORLD)
    :return: Domain object
    """
    return __ReadMesh_driver(filename, integrationOrder, reducedIntegrationOrder,
                             optimize, diracPoints, diracTags, comm)

def ReadGmsh(fileName, numDim, integrationOrder=-1, reducedIntegrationOrder=-1, optimize=True,
             useMacroElements=False, diracPoints=[], diracTags=[], comm=None):
    """
    Read a gmsh mesh file.

    :param fileName: Path to gmsh file
    :param numDim: Number of spatial dimensions
    :param integrationOrder: Order of quadrature scheme (-1 for automatic)
    :param reducedIntegrationOrder: Order of reduced quadrature scheme (-1 for automatic)
    :param optimize: Enable optimization of node labels
    :param useMacroElements: Enable usage of macro elements instead of second order elements
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, defaults to MPI_COMM_WORLD)
    :return: Domain object
    """
    return __ReadGmsh_driver(fileName, numDim, integrationOrder, reducedIntegrationOrder,
                             optimize, useMacroElements, diracPoints, diracTags, comm)


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
    if 'comm' in kwargs:
        mpi=kwargs['comm']
        if mpi=='None':
            args=[n0, n1, order, l0, l1, periodic0, periodic1, integrationOrder, 
                reducedIntegrationOrder, faceon, useFullElementOrder, optimize, points, tags];
            return __Rectangle_driver(args)
        else:
            args=[n0, n1, order, l0, l1, periodic0, periodic1, integrationOrder, 
                reducedIntegrationOrder, faceon, useFullElementOrder, optimize, points, tags, mpi];
            return __Rectangle_driver_MPI(args)    
    else:
        args=[n0, n1, order, l0, l1, periodic0, periodic1, integrationOrder, 
            reducedIntegrationOrder, faceon, useFullElementOrder, optimize, points, tags];
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
    if 'comm' in kwargs:
        mpi=kwargs['comm']
        if mpi=='None':
            args=[n0, n1, n2, order, l0, l1, l2, periodic0,  periodic1, periodic2,
                integrationOrder, reducedIntegrationOrder, faceon, useFullElementOrder,
                optimize, points, tags]
            return __Brick_driver(args)
        else:
            args=[n0, n1, n2, order, l0, l1, l2, periodic0,  periodic1, periodic2,
                integrationOrder, reducedIntegrationOrder, faceon, useFullElementOrder,
                optimize, points, tags, mpi]
            return __Brick_driver_MPI(args)    
    else:
        args=[n0, n1, n2, order, l0, l1, l2, periodic0,  periodic1, periodic2,
            integrationOrder, reducedIntegrationOrder, faceon, useFullElementOrder,
            optimize, points, tags]
        return __Brick_driver(args)
    # args=[n0, n1, n2, order, l0, l1, l2, periodic0,  periodic1, periodic2,
    # integrationOrder, reducedIntegrationOrder, faceon, useFullElementOrder,
    # optimize, points, tags];
    # return __Brick_driver(args)

Brick.__doc__=__Brick_driver.__doc__
