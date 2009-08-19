
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

## :file symbols.py

"""
some tools supporting the usage of symbols.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from util import Symbol

def ScalarSymbol(dim=None):
      """
      Returns a rank 0 L{Symbol}.

      :param dim: spatial dimension or an object that has the C{getDim} method
                  defining the spatial dimension. If dim=C{None}, the spatial
                  diminsion of the returned L{Symbol} is undefined.
      :type dim: C{None}, C{int} or any object with a C{getDim} method
      :return: a L{Symbol} of rank 0
      :rtype: L{Symbol}
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(),dim=d,args=[])


def VectorSymbol(dim=3):
      """
      Returns a vector L{Symbol} of rank 1 and spatial dimension C{dim}.

      :param dim: spatial dimension or an object that has the C{getDim} method
                  defining the spatial dimension
      :type dim: C{int} or any object with a C{getDim} method
      :return: a L{Symbol} of shape (C{dim},)
      :rtype: L{Symbol}
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,),dim=d,args=[])

def TensorSymbol(dim=3):
      """
      Returns a tensor L{Symbol} of rank 2 and spatial dimension C{dim}.

      :param dim: spatial dimension or an object that has the C{getDim} method
                  defining the spatial dimension
      :type dim: C{int} or any object with a C{getDim} method
      :return: a L{Symbol} of shape (C{dim},C{dim})
      :rtype: L{Symbol}
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,d),dim=d,args=[])

def Tensor3Symbol(dim=3):
      """
      Returns a tensor L{Symbol} of rank 3 and spatial dimension C{dim}.

      :param dim: spatial dimension or an object that has the C{getDim} method
                  defining the spatial dimension
      :type dim: C{int} or any object with a C{getDim} method
      :return: a L{Symbol} of shape (C{dim},C{dim},C{dim})
      :rtype: L{Symbol}
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,d,d),dim=d,args=[])

def Tensor4Symbol(dim=3):
      """
      Returns a tensor L{Symbol} of rank 4 and spatial dimension C{dim}.

      :param dim: spatial dimension or an object that has the C{getDim} method
                  defining the spatial dimension
      :type dim: C{int} or any object with a C{getDim} method
      :return: a L{Symbol} of shape (C{dim},C{dim},C{dim},C{dim})
      :rtype: L{Symbol}
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,d,d,d),dim=d,args=[])
#
# $Log:$
#
# vim: expandtab shiftwidth=4:
