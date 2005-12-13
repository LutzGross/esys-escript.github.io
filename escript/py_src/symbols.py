# $Id$

#
#      COPYRIGHT ACcESS 2004 -  All Rights Reserved
#
#   This software is the property of ACcESS.  No part of this code
#   may be copied in any form or by any means without the expressed written
#   consent of ACcESS.  Copying, use or modification of this software
#   by any unauthorised person is illegal unless that
#   person has a software license agreement with ACcESS.
#

## @file symbols.py

"""
some tools supporting the usage of symbols.

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""
                                                                                                                                                                                                     
__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from util import Symbol

def ScalarSymbol(dim=None):
      """
      returns a rank 0 L{Symbol}.

      @param dim: spatial dimension or an object that has the C{getDim} method defining the spatial dimension. If dim=C{None}, the spatial diminsion of the returned L{Symbol} is undefined.
      @type dim: C{None}, C{int} or any object with a C{getDim} method
      @return: a L{Symbol} of rank 0.
      @rtype: L{Symbol} 
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(),dim=d,args=[])


def VectorSymbol(dim=3):
      """
      returns a vector L{Symbol} of rank 1 and spatial dimension C{dim}  

      @param dim: spatial dimension or an object that has the C{getDim} method defining the spatial dimension.
      @type dim: C{int} or any object with a C{getDim} method
      @return: a L{Symbol} of shape (C{dim},) 
      @rtype: L{Symbol} 
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,),dim=d,args=[])

def TensorSymbol(dim=3):
      """
      returns a tensor L{Symbol} of rank 2 and spatial dimension C{dim}

      @param dim: spatial dimension or an object that has the C{getDim} method defining the spatial dimension.
      @type dim: C{int} or any object with a C{getDim} method
      @return: a L{Symbol} of shape (C{dim},C{dim}) 
      @rtype: L{Symbol} 
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:
           d=dim
      return Symbol(shape=(d,d),dim=d,args=[])

def Tensor3Symbol(dim=3):
      """
      returns a tensor L{Symbol} of rank 3 and spatial dimension C{dim}

      @param dim: spatial dimension or an object that has the C{getDim} method defining the spatial dimension.
      @type dim: C{int} or any object with a C{getDim} method
      @return: a L{Symbol} of shape (C{dim},C{dim},C{dim}) 
      @rtype: L{Symbol} 
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      return Symbol(shape=(d,d,d),dim=d,args=[])

def Tensor4Symbol(dim=3):
      """
      returns a tensor L{Symbol} of rank 4 and spatial dimension C{dim}

      @param dim: spatial dimension or an object that has the C{getDim} method defining the spatial dimension.
      @type dim: C{int} or any object with a C{getDim} method
      @param name: name of the symbol
      @type name: C{str}
      @return: a L{Symbol} of shape (C{dim},C{dim},C{dim},C{dim}) 
      @rtype: L{Symbol} 
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
