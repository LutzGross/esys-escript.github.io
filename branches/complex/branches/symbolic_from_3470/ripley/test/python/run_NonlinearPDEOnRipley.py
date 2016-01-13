
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
"""

import sys
import os
import unittest
from esys.escript import *
from esys.ripley import Rectangle,Brick
from esys.escript.linearPDEs import IllegalCoefficient, IllegalCoefficient

# from xx import Test_NonlinearPDE
DEBUG=False

NE=10 # number of elements in each spatial direction

class Test_NonlinearPDE(unittest.TestCase):
   def test_NonLinearPDE_Unknown1_X_Component_0_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (DIM,), "wrong shape of coefficient X.")
       f[:]=g*u
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=u*numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref=g
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown1_X_reduced_Component_0_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (DIM,), "wrong shape of coefficient X_reduced.")
       f[:]=g*u
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=u*numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref=g
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown1_Y_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), () ,"wrong shape of coefficient Y.")
       f=length(g)**2*u
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=length(g)**2
       C_ref=2*g*u
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown1_Y_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), () ,"wrong shape of coefficient Y_reduced.")
       f=length(g)**2*u
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=length(g)**2
       C_ref=2*g*u
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown1_y_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), () ,"wrong shape of coefficient y.")
       f=u
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (), dim=DIM)*0+1
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown1_y_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), () ,"wrong shape of coefficient y_reduced.")
       f=u
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (), dim=DIM)*0+1
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown1_y_contact_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), () ,"wrong shape of coefficient y_contact.")
       f=u
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (), dim=DIM)*0+1
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown1_y_contact_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), () ,"wrong shape of coefficient y_contact_reduced.")
       f=u
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (), dim=DIM)*0+1
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown2_X_Component_0_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown2_X_Component_0_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[0,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown2_X_Component_1_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[1,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[0,:]*u[1]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown2_X_Component_1_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[1,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,1]=2*grad(u)[1,:]*u[1]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown2_X_reduced_Component_0_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown2_X_reduced_Component_0_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[0,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown2_X_reduced_Component_1_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[1,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[0,:]*u[1]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown2_X_reduced_Component_1_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (2,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (2,DIM,2,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (2,DIM,2), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[1,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,1]=2*grad(u)[1,:]*u[1]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown2_Y_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient Y.")
       f[0]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (2,2), dim=DIM)
       C_ref=Symbol('C_ref', (2,2,DIM), dim=DIM)
       D_ref[0,0]=2*length(g)**2*u[0]
       C_ref[0,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[0,1]=2*length(g)**2*u[1]
       C_ref[0,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown2_Y_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient Y.")
       f[1]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (2,2), dim=DIM)
       C_ref=Symbol('C_ref', (2,2,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[1,0]=2*length(g)**2*u[0]
       C_ref[1,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[1,1]=2*length(g)**2*u[1]
       C_ref[1,1,:]=2*length(u)**2*grad(u)[1,:]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown2_Y_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient Y_reduced.")
       f[0]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (2,2), dim=DIM)
       C_ref=Symbol('C_ref', (2,2,DIM), dim=DIM)
       D_ref[0,0]=2*length(g)**2*u[0]
       C_ref[0,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[0,1]=2*length(g)**2*u[1]
       C_ref[0,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown2_Y_reduced_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient Y_reduced.")
       f[1]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (2,2), dim=DIM)
       C_ref=Symbol('C_ref', (2,2,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[1,0]=2*length(g)**2*u[0]
       C_ref[1,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[1,1]=2*length(g)**2*u[1]
       C_ref[1,1,:]=2*length(u)**2*grad(u)[1,:]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown2_y_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y.")
       f[0]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[1,0]=0
       D_ref[1,1]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown2_y_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y.")
       f[1]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown2_y_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y_reduced.")
       f[0]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[1,0]=0
       D_ref[1,1]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown2_y_reduced_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y_reduced.")
       f[1]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown2_y_contact_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y_contact.")
       f[0]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[1,0]=0
       D_ref[1,1]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown2_y_contact_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y_contact.")
       f[1]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown2_y_contact_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y_contact_reduced.")
       f[0]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[1,0]=0
       D_ref[1,1]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown2_y_contact_reduced_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (2,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (2,), "wrong shape of coefficient y_contact_reduced.")
       f[1]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (2, 2), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown5_X_Component_0_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_0_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[0,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_0_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[2,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[0,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_0_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[3,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[0,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_0_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[0,:]=g[4,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[0,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_1_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_1_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[1,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_1_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[2,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[1,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_1_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[3,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[1,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_1_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[1,:]=g[4,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[1,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_2_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[2,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_2_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[2,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[2,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_2_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[2,:]=g[2,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[2,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_2_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[2,:]=g[3,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[2,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_2_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[2,:]=g[4,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[2,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_3_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[3,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_3_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[3,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[3,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_3_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[3,:]=g[2,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[3,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_3_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[3,:]=g[3,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[3,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_3_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[3,:]=g[4,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[3,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_4_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[4,:]=g[0,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[0,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_4_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[4,:]=g[1,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[4,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[1,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_4_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[4,:]=g[2,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[4,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[2,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_4_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[4,:]=g[3,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[4,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[3,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_Component_4_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X.")
       f[4,:]=g[4,:]*length(u)**2
       pde.setValue(X=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X")
       self.assertEqual(X, f, "wrong coefficient X.")
       A=pde.getCoefficient("A")
       B=pde.getCoefficient("B")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[4,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,4]=2*grad(u)[4,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_0_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_0_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[0,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_0_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[2,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[0,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_0_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[3,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[0,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[0,:,4,:] = 0
       B_ref[0,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_0_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[0,:]=g[4,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:] = 0
       B_ref[0,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[0,:,1,:] = 0
       B_ref[0,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[0,:,2,:] = 0
       B_ref[0,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[0,:,3,:] = 0
       B_ref[0,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[0,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[0,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_1_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_1_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[1,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_1_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[2,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[1,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_1_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[3,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[1,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[1,:,4,:] = 0
       B_ref[1,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_1_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[1,:]=g[4,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:] = 0
       B_ref[1,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[1,:,1,:] = 0
       B_ref[1,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[1,:,2,:] = 0
       B_ref[1,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[1,:,3,:] = 0
       B_ref[1,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[1,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[1,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_2_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[2,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_2_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[2,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[2,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_2_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[2,:]=g[2,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[2,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_2_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[2,:]=g[3,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[2,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[2,:,4,:] = 0
       B_ref[2,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_2_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[2,:]=g[4,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:] = 0
       B_ref[2,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[2,:,1,:] = 0
       B_ref[2,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[2,:,2,:] = 0
       B_ref[2,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[2,:,3,:] = 0
       B_ref[2,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[2,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[2,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_3_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[3,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[0,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_3_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[3,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[3,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[1,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_3_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[3,:]=g[2,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[3,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[2,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_3_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[3,:]=g[3,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[3,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[3,:,4,:] = 0
       B_ref[3,:,4]=2*grad(u)[3,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_3_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[3,:]=g[4,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:] = 0
       B_ref[3,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[3,:,1,:] = 0
       B_ref[3,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[3,:,2,:] = 0
       B_ref[3,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[3,:,3,:] = 0
       B_ref[3,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[3,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[3,:,4]=2*grad(u)[4,:]*u[4]
       A_ref[4,:,0,:]=0
       B_ref[4,:, 0]=0
       A_ref[4,:,1,:]=0
       B_ref[4,:, 1]=0
       A_ref[4,:,2,:]=0
       B_ref[4,:, 2]=0
       A_ref[4,:,3,:]=0
       B_ref[4,:, 3]=0
       A_ref[4,:,4,:]=0
       B_ref[4,:, 4]=0
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_4_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[4,:]=g[0,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,0]=2*grad(u)[0,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[0,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[0,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[0,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[0,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_4_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[4,:]=g[1,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[1,:]*u[0]
       A_ref[4,:,1,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,1]=2*grad(u)[1,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[1,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[1,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[1,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_4_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[4,:]=g[2,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[2,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[2,:]*u[1]
       A_ref[4,:,2,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,2]=2*grad(u)[2,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[2,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[2,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_4_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[4,:]=g[3,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[3,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[3,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[3,:]*u[2]
       A_ref[4,:,3,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,3]=2*grad(u)[3,:]*u[3]
       A_ref[4,:,4,:] = 0
       B_ref[4,:,4]=2*grad(u)[3,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_X_reduced_Component_4_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("X_reduced")
       self.assertEqual(f.getShape(), (5,DIM), "wrong shape of coefficient X_reduced.")
       f[4,:]=g[4,:]*length(u)**2
       pde.setValue(X_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       X=pde.getCoefficient("X_reduced")
       self.assertEqual(X, f, "wrong coefficient X_reduced.")
       A=pde.getCoefficient("A_reduced")
       B=pde.getCoefficient("B_reduced")
       A_ref=Symbol('A_ref', (5,DIM,5,DIM), dim=DIM)
       B_ref=Symbol('B_ref', (5,DIM,5), dim=DIM)
       A_ref[0,:,0,:]=0
       B_ref[0,:, 0]=0
       A_ref[0,:,1,:]=0
       B_ref[0,:, 1]=0
       A_ref[0,:,2,:]=0
       B_ref[0,:, 2]=0
       A_ref[0,:,3,:]=0
       B_ref[0,:, 3]=0
       A_ref[0,:,4,:]=0
       B_ref[0,:, 4]=0
       A_ref[1,:,0,:]=0
       B_ref[1,:, 0]=0
       A_ref[1,:,1,:]=0
       B_ref[1,:, 1]=0
       A_ref[1,:,2,:]=0
       B_ref[1,:, 2]=0
       A_ref[1,:,3,:]=0
       B_ref[1,:, 3]=0
       A_ref[1,:,4,:]=0
       B_ref[1,:, 4]=0
       A_ref[2,:,0,:]=0
       B_ref[2,:, 0]=0
       A_ref[2,:,1,:]=0
       B_ref[2,:, 1]=0
       A_ref[2,:,2,:]=0
       B_ref[2,:, 2]=0
       A_ref[2,:,3,:]=0
       B_ref[2,:, 3]=0
       A_ref[2,:,4,:]=0
       B_ref[2,:, 4]=0
       A_ref[3,:,0,:]=0
       B_ref[3,:, 0]=0
       A_ref[3,:,1,:]=0
       B_ref[3,:, 1]=0
       A_ref[3,:,2,:]=0
       B_ref[3,:, 2]=0
       A_ref[3,:,3,:]=0
       B_ref[3,:, 3]=0
       A_ref[3,:,4,:]=0
       B_ref[3,:, 4]=0
       A_ref[4,:,0,:] = 0
       B_ref[4,:,0]=2*grad(u)[4,:]*u[0]
       A_ref[4,:,1,:] = 0
       B_ref[4,:,1]=2*grad(u)[4,:]*u[1]
       A_ref[4,:,2,:] = 0
       B_ref[4,:,2]=2*grad(u)[4,:]*u[2]
       A_ref[4,:,3,:] = 0
       B_ref[4,:,3]=2*grad(u)[4,:]*u[3]
       A_ref[4,:,4,:] = length(u)**2 *numpy.array(kronecker(DIM), dtype=numpy.int32)
       B_ref[4,:,4]=2*grad(u)[4,:]*u[4]
       self.assertEqual(B_ref.simplify(), B.simplify(), "wrong coefficient B_reduced.")
       self.assertEqual(A_ref.simplify(), A.simplify(), "wrong coefficient A_reduced.")
   def test_NonLinearPDE_Unknown5_Y_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y.")
       f[0]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=2*length(g)**2*u[0]
       C_ref[0,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[0,1]=2*length(g)**2*u[1]
       C_ref[0,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[0,2]=2*length(g)**2*u[2]
       C_ref[0,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[0,3]=2*length(g)**2*u[3]
       C_ref[0,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[0,4]=2*length(g)**2*u[4]
       C_ref[0,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown5_Y_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y.")
       f[1]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=2*length(g)**2*u[0]
       C_ref[1,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[1,1]=2*length(g)**2*u[1]
       C_ref[1,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[1,2]=2*length(g)**2*u[2]
       C_ref[1,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[1,3]=2*length(g)**2*u[3]
       C_ref[1,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[1,4]=2*length(g)**2*u[4]
       C_ref[1,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown5_Y_Component_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y.")
       f[2]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=2*length(g)**2*u[0]
       C_ref[2,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[2,1]=2*length(g)**2*u[1]
       C_ref[2,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[2,2]=2*length(g)**2*u[2]
       C_ref[2,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[2,3]=2*length(g)**2*u[3]
       C_ref[2,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[2,4]=2*length(g)**2*u[4]
       C_ref[2,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown5_Y_Component_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y.")
       f[3]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=2*length(g)**2*u[0]
       C_ref[3,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[3,1]=2*length(g)**2*u[1]
       C_ref[3,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[3,2]=2*length(g)**2*u[2]
       C_ref[3,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[3,3]=2*length(g)**2*u[3]
       C_ref[3,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[3,4]=2*length(g)**2*u[4]
       C_ref[3,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown5_Y_Component_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y.")
       f[4]=length(g)**2*length(u)**2
       pde.setValue(Y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y")
       self.assertEqual(Y, f, "wrong coefficient Y.")
       D=pde.getCoefficient("D")
       C=pde.getCoefficient("C")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=2*length(g)**2*u[0]
       C_ref[4,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[4,1]=2*length(g)**2*u[1]
       C_ref[4,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[4,2]=2*length(g)**2*u[2]
       C_ref[4,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[4,3]=2*length(g)**2*u[3]
       C_ref[4,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[4,4]=2*length(g)**2*u[4]
       C_ref[4,4,:]=2*length(u)**2*grad(u)[4,:]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C.")
   def test_NonLinearPDE_Unknown5_Y_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y_reduced.")
       f[0]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=2*length(g)**2*u[0]
       C_ref[0,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[0,1]=2*length(g)**2*u[1]
       C_ref[0,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[0,2]=2*length(g)**2*u[2]
       C_ref[0,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[0,3]=2*length(g)**2*u[3]
       C_ref[0,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[0,4]=2*length(g)**2*u[4]
       C_ref[0,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown5_Y_reduced_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y_reduced.")
       f[1]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=2*length(g)**2*u[0]
       C_ref[1,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[1,1]=2*length(g)**2*u[1]
       C_ref[1,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[1,2]=2*length(g)**2*u[2]
       C_ref[1,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[1,3]=2*length(g)**2*u[3]
       C_ref[1,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[1,4]=2*length(g)**2*u[4]
       C_ref[1,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown5_Y_reduced_Component_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y_reduced.")
       f[2]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=2*length(g)**2*u[0]
       C_ref[2,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[2,1]=2*length(g)**2*u[1]
       C_ref[2,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[2,2]=2*length(g)**2*u[2]
       C_ref[2,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[2,3]=2*length(g)**2*u[3]
       C_ref[2,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[2,4]=2*length(g)**2*u[4]
       C_ref[2,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown5_Y_reduced_Component_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y_reduced.")
       f[3]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=2*length(g)**2*u[0]
       C_ref[3,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[3,1]=2*length(g)**2*u[1]
       C_ref[3,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[3,2]=2*length(g)**2*u[2]
       C_ref[3,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[3,3]=2*length(g)**2*u[3]
       C_ref[3,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[3,4]=2*length(g)**2*u[4]
       C_ref[3,4,:]=2*length(u)**2*grad(u)[4,:]
       D_ref[4,0]=0
       C_ref[4,0,:]=0
       D_ref[4,1]=0
       C_ref[4,1,:]=0
       D_ref[4,2]=0
       C_ref[4,2,:]=0
       D_ref[4,3]=0
       C_ref[4,3,:]=0
       D_ref[4,4]=0
       C_ref[4,4,:]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown5_Y_reduced_Component_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("Y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient Y_reduced.")
       f[4]=length(g)**2*length(u)**2
       pde.setValue(Y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("Y_reduced")
       self.assertEqual(Y, f, "wrong coefficient Y_reduced.")
       D=pde.getCoefficient("D_reduced")
       C=pde.getCoefficient("C_reduced")
       D_ref=Symbol('D_ref', (5,5), dim=DIM)
       C_ref=Symbol('C_ref', (5,5,DIM), dim=DIM)
       D_ref[0,0]=0
       C_ref[0,0,:]=0
       D_ref[0,1]=0
       C_ref[0,1,:]=0
       D_ref[0,2]=0
       C_ref[0,2,:]=0
       D_ref[0,3]=0
       C_ref[0,3,:]=0
       D_ref[0,4]=0
       C_ref[0,4,:]=0
       D_ref[1,0]=0
       C_ref[1,0,:]=0
       D_ref[1,1]=0
       C_ref[1,1,:]=0
       D_ref[1,2]=0
       C_ref[1,2,:]=0
       D_ref[1,3]=0
       C_ref[1,3,:]=0
       D_ref[1,4]=0
       C_ref[1,4,:]=0
       D_ref[2,0]=0
       C_ref[2,0,:]=0
       D_ref[2,1]=0
       C_ref[2,1,:]=0
       D_ref[2,2]=0
       C_ref[2,2,:]=0
       D_ref[2,3]=0
       C_ref[2,3,:]=0
       D_ref[2,4]=0
       C_ref[2,4,:]=0
       D_ref[3,0]=0
       C_ref[3,0,:]=0
       D_ref[3,1]=0
       C_ref[3,1,:]=0
       D_ref[3,2]=0
       C_ref[3,2,:]=0
       D_ref[3,3]=0
       C_ref[3,3,:]=0
       D_ref[3,4]=0
       C_ref[3,4,:]=0
       D_ref[4,0]=2*length(g)**2*u[0]
       C_ref[4,0,:]=2*length(u)**2*grad(u)[0,:]
       D_ref[4,1]=2*length(g)**2*u[1]
       C_ref[4,1,:]=2*length(u)**2*grad(u)[1,:]
       D_ref[4,2]=2*length(g)**2*u[2]
       C_ref[4,2,:]=2*length(u)**2*grad(u)[2,:]
       D_ref[4,3]=2*length(g)**2*u[3]
       C_ref[4,3,:]=2*length(u)**2*grad(u)[3,:]
       D_ref[4,4]=2*length(g)**2*u[4]
       C_ref[4,4,:]=2*length(u)**2*grad(u)[4,:]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
       self.assertEqual(C_ref.simplify(), C.simplify(), "wrong coefficient C_reduced.")
   def test_NonLinearPDE_Unknown5_y_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y.")
       f[0]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[0,2]=2*u[2]
       D_ref[0,3]=2*u[3]
       D_ref[0,4]=2*u[4]
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown5_y_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y.")
       f[1]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       D_ref[1,2]=2*u[2]
       D_ref[1,3]=2*u[3]
       D_ref[1,4]=2*u[4]
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown5_y_Component_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y.")
       f[2]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=2*u[0]
       D_ref[2,1]=2*u[1]
       D_ref[2,2]=2*u[2]
       D_ref[2,3]=2*u[3]
       D_ref[2,4]=2*u[4]
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown5_y_Component_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y.")
       f[3]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=2*u[0]
       D_ref[3,1]=2*u[1]
       D_ref[3,2]=2*u[2]
       D_ref[3,3]=2*u[3]
       D_ref[3,4]=2*u[4]
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown5_y_Component_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y.")
       f[4]=length(u)**2
       pde.setValue(y=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y")
       self.assertEqual(Y, f, "wrong coefficient y.")
       D=pde.getCoefficient("d")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=2*u[0]
       D_ref[4,1]=2*u[1]
       D_ref[4,2]=2*u[2]
       D_ref[4,3]=2*u[3]
       D_ref[4,4]=2*u[4]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D.")
   def test_NonLinearPDE_Unknown5_y_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_reduced.")
       f[0]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[0,2]=2*u[2]
       D_ref[0,3]=2*u[3]
       D_ref[0,4]=2*u[4]
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown5_y_reduced_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_reduced.")
       f[1]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       D_ref[1,2]=2*u[2]
       D_ref[1,3]=2*u[3]
       D_ref[1,4]=2*u[4]
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown5_y_reduced_Component_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_reduced.")
       f[2]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=2*u[0]
       D_ref[2,1]=2*u[1]
       D_ref[2,2]=2*u[2]
       D_ref[2,3]=2*u[3]
       D_ref[2,4]=2*u[4]
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown5_y_reduced_Component_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_reduced.")
       f[3]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=2*u[0]
       D_ref[3,1]=2*u[1]
       D_ref[3,2]=2*u[2]
       D_ref[3,3]=2*u[3]
       D_ref[3,4]=2*u[4]
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown5_y_reduced_Component_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_reduced.")
       f[4]=length(u)**2
       pde.setValue(y_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_reduced.")
       D=pde.getCoefficient("d_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=2*u[0]
       D_ref[4,1]=2*u[1]
       D_ref[4,2]=2*u[2]
       D_ref[4,3]=2*u[3]
       D_ref[4,4]=2*u[4]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_reduced.")
   def test_NonLinearPDE_Unknown5_y_contact_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact.")
       f[0]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[0,2]=2*u[2]
       D_ref[0,3]=2*u[3]
       D_ref[0,4]=2*u[4]
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown5_y_contact_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact.")
       f[1]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       D_ref[1,2]=2*u[2]
       D_ref[1,3]=2*u[3]
       D_ref[1,4]=2*u[4]
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown5_y_contact_Component_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact.")
       f[2]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=2*u[0]
       D_ref[2,1]=2*u[1]
       D_ref[2,2]=2*u[2]
       D_ref[2,3]=2*u[3]
       D_ref[2,4]=2*u[4]
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown5_y_contact_Component_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact.")
       f[3]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=2*u[0]
       D_ref[3,1]=2*u[1]
       D_ref[3,2]=2*u[2]
       D_ref[3,3]=2*u[3]
       D_ref[3,4]=2*u[4]
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown5_y_contact_Component_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact.")
       f[4]=length(u)**2
       pde.setValue(y_contact=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact_reduced")
       Y=pde.getCoefficient("y_contact")
       self.assertEqual(Y, f, "wrong coefficient y_contact.")
       D=pde.getCoefficient("d_contact")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=2*u[0]
       D_ref[4,1]=2*u[1]
       D_ref[4,2]=2*u[2]
       D_ref[4,3]=2*u[3]
       D_ref[4,4]=2*u[4]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact.")
   def test_NonLinearPDE_Unknown5_y_contact_reduced_Component_0(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact_reduced.")
       f[0]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=2*u[0]
       D_ref[0,1]=2*u[1]
       D_ref[0,2]=2*u[2]
       D_ref[0,3]=2*u[3]
       D_ref[0,4]=2*u[4]
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown5_y_contact_reduced_Component_1(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact_reduced.")
       f[1]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=2*u[0]
       D_ref[1,1]=2*u[1]
       D_ref[1,2]=2*u[2]
       D_ref[1,3]=2*u[3]
       D_ref[1,4]=2*u[4]
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown5_y_contact_reduced_Component_2(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact_reduced.")
       f[2]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=2*u[0]
       D_ref[2,1]=2*u[1]
       D_ref[2,2]=2*u[2]
       D_ref[2,3]=2*u[3]
       D_ref[2,4]=2*u[4]
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown5_y_contact_reduced_Component_3(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact_reduced.")
       f[3]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=2*u[0]
       D_ref[3,1]=2*u[1]
       D_ref[3,2]=2*u[2]
       D_ref[3,3]=2*u[3]
       D_ref[3,4]=2*u[4]
       D_ref[4,0]=0
       D_ref[4,1]=0
       D_ref[4,2]=0
       D_ref[4,3]=0
       D_ref[4,4]=0
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")
   def test_NonLinearPDE_Unknown5_y_contact_reduced_Component_4(self):
       DIM = self.domain.getDim()
       u=Symbol('v', (5,), dim=DIM)
       pde=NonlinearPDE(self.domain, u, debug=DEBUG)
       self.assertEqual(pde.getUnknownSymbol(), u, "solution symbol is wrong.")
       g=grad(u)
       f=pde.createCoefficient("y_contact_reduced")
       self.assertEqual(f.getShape(), (5,), "wrong shape of coefficient y_contact_reduced.")
       f[4]=length(u)**2
       pde.setValue(y_contact_reduced=f)
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_contact")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "X_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "A_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "B_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "Y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "D_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "C_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "y_reduced")
       self.assertRaises(IllegalCoefficient, pde.getCoefficient, "d_reduced")
       Y=pde.getCoefficient("y_contact_reduced")
       self.assertEqual(Y, f, "wrong coefficient y_contact_reduced.")
       D=pde.getCoefficient("d_contact_reduced")
       D_ref=Symbol('D_ref', (5, 5), dim=DIM)
       D_ref[0,0]=0
       D_ref[0,1]=0
       D_ref[0,2]=0
       D_ref[0,3]=0
       D_ref[0,4]=0
       D_ref[1,0]=0
       D_ref[1,1]=0
       D_ref[1,2]=0
       D_ref[1,3]=0
       D_ref[1,4]=0
       D_ref[2,0]=0
       D_ref[2,1]=0
       D_ref[2,2]=0
       D_ref[2,3]=0
       D_ref[2,4]=0
       D_ref[3,0]=0
       D_ref[3,1]=0
       D_ref[3,2]=0
       D_ref[3,3]=0
       D_ref[3,4]=0
       D_ref[4,0]=2*u[0]
       D_ref[4,1]=2*u[1]
       D_ref[4,2]=2*u[2]
       D_ref[4,3]=2*u[3]
       D_ref[4,4]=2*u[4]
       self.assertEqual(D_ref.simplify(), D.simplify(), "wrong coefficient D_contact_reduced.")


class Test_NonlinearPDEOnRipley2D(Test_NonlinearPDE):
   def setUp(self):
        self.domain = Rectangle(NE,NE)
   def tearDown(self):
        del self.domain

class Test_NonlinearPDEOnRipley3D(Test_NonlinearPDE):
   def setUp(self):
        self.domain = Brick(NE,NE,NE)
   def tearDown(self):
        del self.domain
        
if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_NonlinearPDEOnRipley2D))
   suite.addTest(unittest.makeSuite(Test_NonlinearPDEOnRipley3D))
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)        
