# $Id$

"""
Test suite for the symbols.py module

@remark:

@var __author__: name of author
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from esys.escript import *
from esys.escript.symbols import *
import unittest

class Test_symbols(unittest.TestCase):

   def test_Symbol_Scalar_dNone(self):
      s=ScalarSymbol(dim=None)
      self.failUnlessEqual(s.getRank(),0,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(),"wrong shape.")
      self.failUnlessEqual(s.getDim(),None,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Scalar_dd(self):
      s=ScalarSymbol(dim=self.functionspace)
      d=self.functionspace.getDim()
      self.failUnlessEqual(s.getRank(),0,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(),"wrong shape.")
      self.failUnlessEqual(s.getDim(),d,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Scalar_d1(self):
      s=ScalarSymbol(dim=1)
      self.failUnlessEqual(s.getRank(),0,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(),"wrong shape.")
      self.failUnlessEqual(s.getDim(),1,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Scalar_d2(self):
      s=ScalarSymbol(dim=2)
      self.failUnlessEqual(s.getRank(),0,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(),"wrong shape.")
      self.failUnlessEqual(s.getDim(),2,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Scalar_d3(self):
      s=ScalarSymbol(dim=3)
      self.failUnlessEqual(s.getRank(),0,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(),"wrong shape.")
      self.failUnlessEqual(s.getDim(),3,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Vector_dd(self):
      s=VectorSymbol(dim=self.functionspace)
      d=self.functionspace.getDim()
      self.failUnlessEqual(s.getRank(),1,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(d,),"wrong shape.")
      self.failUnlessEqual(s.getDim(),d,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Vector_d1(self):
      s=VectorSymbol(dim=1)
      self.failUnlessEqual(s.getRank(),1,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(1,),"wrong shape.")
      self.failUnlessEqual(s.getDim(),1,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Vector_d2(self):
      s=VectorSymbol(dim=2)
      self.failUnlessEqual(s.getRank(),1,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(2,),"wrong shape.")
      self.failUnlessEqual(s.getDim(),2,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Vector_d3(self):
      s=VectorSymbol(dim=3)
      self.failUnlessEqual(s.getRank(),1,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(3,),"wrong shape.")
      self.failUnlessEqual(s.getDim(),3,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor_dd(self):
      s=TensorSymbol(dim=self.functionspace)
      d=self.functionspace.getDim()
      self.failUnlessEqual(s.getRank(),2,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(d,d),"wrong shape.")
      self.failUnlessEqual(s.getDim(),d,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor_d1(self):
      s=TensorSymbol(dim=1)
      self.failUnlessEqual(s.getRank(),2,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(1, 1),"wrong shape.")
      self.failUnlessEqual(s.getDim(),1,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor_d2(self):
      s=TensorSymbol(dim=2)
      self.failUnlessEqual(s.getRank(),2,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(2, 2),"wrong shape.")
      self.failUnlessEqual(s.getDim(),2,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor_d3(self):
      s=TensorSymbol(dim=3)
      self.failUnlessEqual(s.getRank(),2,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(3, 3),"wrong shape.")
      self.failUnlessEqual(s.getDim(),3,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor3_dd(self):
      s=Tensor3Symbol(dim=self.functionspace)
      d=self.functionspace.getDim()
      self.failUnlessEqual(s.getRank(),3,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(d,d,d),"wrong shape.")
      self.failUnlessEqual(s.getDim(),d,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor3_d1(self):
      s=Tensor3Symbol(dim=1)
      self.failUnlessEqual(s.getRank(),3,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(1, 1, 1),"wrong shape.")
      self.failUnlessEqual(s.getDim(),1,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor3_d2(self):
      s=Tensor3Symbol(dim=2)
      self.failUnlessEqual(s.getRank(),3,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(2, 2, 2),"wrong shape.")
      self.failUnlessEqual(s.getDim(),2,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor3_d3(self):
      s=Tensor3Symbol(dim=3)
      self.failUnlessEqual(s.getRank(),3,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(3, 3, 3),"wrong shape.")
      self.failUnlessEqual(s.getDim(),3,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor4_dd(self):
      s=Tensor4Symbol(dim=self.functionspace)
      d=self.functionspace.getDim()
      self.failUnlessEqual(s.getRank(),4,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(d,d,d,d),"wrong shape.")
      self.failUnlessEqual(s.getDim(),d,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor4_d1(self):
      s=Tensor4Symbol(dim=1)
      self.failUnlessEqual(s.getRank(),4,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(1, 1, 1, 1),"wrong shape.")
      self.failUnlessEqual(s.getDim(),1,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor4_d2(self):
      s=Tensor4Symbol(dim=2)
      self.failUnlessEqual(s.getRank(),4,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(2, 2, 2, 2),"wrong shape.")
      self.failUnlessEqual(s.getDim(),2,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Symbol_Tensor4_d3(self):
      s=Tensor4Symbol(dim=3)
      self.failUnlessEqual(s.getRank(),4,"wrong rank.")
      self.failUnlessEqual(s.getShape(),(3, 3, 3, 3),"wrong shape.")
      self.failUnlessEqual(s.getDim(),3,"wrong spatial dimension.")
      self.failUnlessEqual(s.getArgument(),[],"wrong arguments.")
