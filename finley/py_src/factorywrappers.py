
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"


from .finleycpp import __Brick_driver, __Rectangle_driver, __ReadMesh_driver, __ReadGmsh_driver


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
               reducedIntegrationOrder=-1, useElementsOnFace=None, useFullElementOrder=False, optimize=False,
               diracPoints=[], diracTags=[], comm=None):
    """
    Creates a rectangular mesh with n0 x n1 elements over the rectangle [0,l0] x [0,l1].

    :param n0: number of elements in direction 0
    :param n1: number of elements in direction 1
    :param order: order of shape functions (1, 2, or -1 for macro elements)
    :param l0: length of side 0
    :param l1: length of side 1
    :param periodic0: periodic boundary conditions in direction 0
    :param periodic1: periodic boundary conditions in direction 1
    :param integrationOrder: order of quadrature scheme (-1 for automatic)
    :param reducedIntegrationOrder: order of reduced quadrature scheme (-1 for automatic)
    :param useElementsOnFace: whether to use elements on face
    :param useFullElementOrder: whether to use Rec9 elements
    :param optimize: enable optimization of node labels
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, defaults to MPI_COMM_WORLD)
    :return: Domain object
    """
    faceon = useElementsOnFace
    if useElementsOnFace is None:  # We want to use True as the default, but only where it makes sense
        if useFullElementOrder or order == -1:
            faceon = False  # Don't use it
        else:
            faceon = True

    return __Rectangle_driver((n0, n1), order, (l0, l1), (periodic0, periodic1),
                              integrationOrder, reducedIntegrationOrder,
                              faceon, useFullElementOrder, optimize,
                              diracPoints, diracTags, comm)

def Brick(n0=1, n1=1, n2=1, order=1, l0=1.0, l1=1.0, l2=1.0, periodic0=False, periodic1=False, periodic2=False,
          integrationOrder=-1, reducedIntegrationOrder=-1, useElementsOnFace=None, useFullElementOrder=False,
          optimize=False, diracPoints=[], diracTags=[], comm=None):
    """
    Creates a rectangular mesh with n0 x n1 x n2 elements over the brick [0,l0] x [0,l1] x [0,l2].

    :param n0: number of elements in direction 0
    :param n1: number of elements in direction 1
    :param n2: number of elements in direction 2
    :param order: order of shape functions (1, 2, or -1 for macro elements)
    :param l0: length of side 0
    :param l1: length of side 1
    :param l2: length of side 2
    :param periodic0: periodic boundary conditions in direction 0
    :param periodic1: periodic boundary conditions in direction 1
    :param periodic2: periodic boundary conditions in direction 2
    :param integrationOrder: order of quadrature scheme (-1 for automatic)
    :param reducedIntegrationOrder: order of reduced quadrature scheme (-1 for automatic)
    :param useElementsOnFace: whether to use elements on face
    :param useFullElementOrder: whether to use Hex27 elements
    :param optimize: enable optimization of node labels
    :param diracPoints: Dirac point coordinates
    :param diracTags: Dirac point tags
    :param comm: MPI communicator (optional, defaults to MPI_COMM_WORLD)
    :return: Domain object
    """
    faceon = useElementsOnFace
    if useElementsOnFace is None:  # We want to use True as the default, but only where it makes sense
        if useFullElementOrder or order == -1:
            faceon = False  # Don't use it
        else:
            faceon = True

    return __Brick_driver((n0, n1, n2), order, (l0, l1, l2), (periodic0, periodic1, periodic2),
                          integrationOrder, reducedIntegrationOrder,
                          faceon, useFullElementOrder, optimize,
                          diracPoints, diracTags, comm)
