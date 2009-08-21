
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
      Returns a rank 0 `Symbol`.

      :param dim: spatial dimension or an object that has the ``getDim`` method
                  defining the spatial dimension. If dim=``None``, the spatial
                  diminsion of the returned `Symbol` is undefined.
      :type dim: ``None``, ``int`` or any object with a ``getDim`` method
      :return: a `Symbol` of rank 0
      :rtype: `Symbol`
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(),dim=d,args=[])


def VectorSymbol(dim=3):
      """
      Returns a vector `Symbol` of rank 1 and spatial dimension ``dim``.

      :param dim: spatial dimension or an object that has the ``getDim`` method
                  defining the spatial dimension
      :type dim: ``int`` or any object with a ``getDim`` method
      :return: a `Symbol` of shape (``dim``,)
      :rtype: `Symbol`
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,),dim=d,args=[])

def TensorSymbol(dim=3):
      """
      Returns a tensor `Symbol` of rank 2 and spatial dimension ``dim``.

      :param dim: spatial dimension or an object that has the ``getDim`` method
                  defining the spatial dimension
      :type dim: ``int`` or any object with a ``getDim`` method
      :return: a `Symbol` of shape (``dim``,``dim``)
      :rtype: `Symbol`
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,d),dim=d,args=[])

def Tensor3Symbol(dim=3):
      """
      Returns a tensor `Symbol` of rank 3 and spatial dimension ``dim``.

      :param dim: spatial dimension or an object that has the ``getDim`` method
                  defining the spatial dimension
      :type dim: ``int`` or any object with a ``getDim`` method
      :return: a `Symbol` of shape (``dim``,``dim``,``dim``)
      :rtype: `Symbol`
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,d,d),dim=d,args=[])

def Tensor4Symbol(dim=3):
      """
      Returns a tensor `Symbol` of rank 4 and spatial dimension ``dim``.

      :param dim: spatial dimension or an object that has the ``getDim`` method
                  defining the spatial dimension
      :type dim: ``int`` or any object with a ``getDim`` method
      :return: a `Symbol` of shape (``dim``,``dim``,``dim``,``dim``)
      :rtype: `Symbol`
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
