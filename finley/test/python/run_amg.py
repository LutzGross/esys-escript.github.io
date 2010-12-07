# -*- coding: utf-8 -*-

########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for AMG

:remark:

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import unittest, sys

from esys.escript import *
from esys.finley import Rectangle,Brick
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import numpy
OPTIMIZE=True
SOLVER_VERBOSE=True

MIN_MATRIX_SIZE=None
MIN_SPARSITY=None
MAX_LEVEL=None
USE_AMG=True or False

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'
FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

# number of elements in the spatial directions
NE_TOTAL=5000

class Test_AMG(unittest.TestCase):

   def test_Poisson(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex+=(i+1)*x[i]
           g_ex[i]=(i+1)

        # create PDE:
        pde=LinearPDE(self.domain,numEquations=1)
        pde.setSymmetryOn()
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(self.domain),y=inner(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_PoissonSqueezedX(self):
        x=self.domain.getX().copy()
        x[0]*=0.5
        self.domain.setX(x)
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex+=(i+1)*x[i]
           g_ex[i]=(i+1)

        # create PDE:
        pde=LinearPDE(self.domain,numEquations=1)
        pde.setSymmetryOn()
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(self.domain),y=inner(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)


   def test_Poisson2(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1) *x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1

        # create PDE:
        pde=LinearPDE(self.domain,numEquations=2)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_Poisson3(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1) *x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1

        # create PDE:
        pde=LinearPDE(self.domain,numEquations=3)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_Poisson4(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1) *x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           u_ex[3]+= 4*(i+1)*x[i]
           g_ex[3,i]=4*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1
           A[3,i,3,i]=1

        # create PDE:
        pde=LinearPDE(self.domain,numEquations=4)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_WeakCoupled2(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1

        Y=Data(0.,(2,),Function(self.domain))
        Y[0]=        u_ex[1]
        Y[1]=u_ex[0]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=2)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=-0.01*(numpy.ones((2,2))-kronecker(2)),
                     Y=-0.01*Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_WeakCoupled3(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1

        Y=Data(0.,(3,),Function(self.domain))
        Y[0]=        u_ex[1]+u_ex[2]
        Y[1]=u_ex[0]+        u_ex[2]
        Y[2]=u_ex[0]+u_ex[1] 
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=3)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=-0.01*(numpy.ones((3,3))-kronecker(3)),
                     Y=-0.01*Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_WeakCoupled4(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           u_ex[3]+= 4*(i+1)*x[i]
           g_ex[3,i]=4*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1
           A[3,i,3,i]=1

        Y=Data(0.,(4,),Function(self.domain))
        Y[0]=        u_ex[1]+u_ex[2]+u_ex[3]
        Y[1]=u_ex[0]+        u_ex[2]+u_ex[3]
        Y[2]=u_ex[0]+u_ex[1]+        u_ex[3]
        Y[3]=u_ex[0]+u_ex[1]+u_ex[2]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=4)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=-0.01*(numpy.ones((4,4))-kronecker(4)),
                     Y=-0.01*Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_StrongCoupled2(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1

        Y=Data(0.,(2,),Function(self.domain))
        Y[0]=        u_ex[1]
        Y[1]=u_ex[0]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=2)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=-20.*(numpy.ones((2,2))-kronecker(2)),
                     Y=-20.*Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

   def test_StrongCoupled3(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1

        Y=Data(0.,(3,),Function(self.domain))
        Y[0]=        u_ex[1]+u_ex[2]
        Y[1]=u_ex[0]+        u_ex[2]
        Y[2]=u_ex[0]+u_ex[1]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=3)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=-20.*(numpy.ones((3,3))-kronecker(3)),
                     Y=-20.*Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)
   def test_StrongCoupled4(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in xrange(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           u_ex[3]+= 4*(i+1)*x[i]
           g_ex[3,i]=4*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1
           A[3,i,3,i]=1

        Y=Data(0.,(4,),Function(self.domain))
        Y[0]=        u_ex[1]+u_ex[2]+u_ex[3]
        Y[1]=u_ex[0]+        u_ex[2]+u_ex[3]
        Y[2]=u_ex[0]+u_ex[1]+        u_ex[3]
        Y[3]=u_ex[0]+u_ex[1]+u_ex[2]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=4)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=-20.*(numpy.ones((4,4))-kronecker(4)),
                     Y=-20.*Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        if (USE_AMG): pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<self.RES_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class Test_AMGOnFinleyHex2DOrder1(Test_AMG):
   RES_TOL=5.e-7
   SOLVER_TOL=1.e-8
   def setUp(self):
        NE=int(float(NE_TOTAL)**(1./2.)+0.5)
        self.domain = Rectangle(NE,NE,1, optimize=OPTIMIZE)
   def tearDown(self):
        del self.domain

class Test_AMGOnFinleyHex3DOrder1(Test_AMG):
   RES_TOL=5.e-7
   SOLVER_TOL=1.e-8
   def setUp(self):
        NE=int(float(NE_TOTAL)**(1./3.)+0.5)
        self.domain = Brick(NE,NE,NE,1, optimize=OPTIMIZE)
   def tearDown(self):
        del self.domain
if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_AMGOnFinleyHex2DOrder1))
   suite.addTest(unittest.makeSuite(Test_AMGOnFinleyHex3DOrder1))
   #suite.addTest(Test_AMGOnFinleyHex2DOrder1("test_PoissonSqueezedX"))

   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)
