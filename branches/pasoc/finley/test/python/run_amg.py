# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
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

import esys.escriptcore.utestselect as unittest, sys
from esys.escriptcore.testing import *
from esys.escript import *
from esys.finley import Rectangle, Brick
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import numpy

OPTIMIZE = True
SOLVER_VERBOSE = False

MIN_MATRIX_SIZE=1
MIN_SPARSITY=1.
MIN_MATRIX_SIZE=None
MIN_SPARSITY=None
MAX_LEVEL=None
HAVE_PASO_AMG = hasFeature('paso') and (getMPISizeWorld() == 1)

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'
FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

# number of elements in the spatial directions
NE_TOTAL=4096
#NE_TOTAL=4

@unittest.skipUnless(HAVE_PASO_AMG, "PASO AMG tests require PASO and #ranks == 1")
class AMG(unittest.TestCase): #subclassing required

   def test_Poisson(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_PoissonWithDirectInterpolation(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setNumPreSweeps(3)
        pde.getSolverOptions().setNumPostSweeps(3)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        pde.getSolverOptions().setAMGInterpolation(SolverOptions.DIRECT_INTERPOLATION)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_PoissonClassic(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        pde.getSolverOptions().setNumPreSweeps(3)
        pde.getSolverOptions().setNumPostSweeps(3)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        pde.getSolverOptions().setAMGInterpolation(SolverOptions.CLASSIC_INTERPOLATION)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_PoissonClassicWithFFCoupling(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        pde.getSolverOptions().setNumPreSweeps(3)
        pde.getSolverOptions().setNumPostSweeps(3)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        pde.getSolverOptions().setAMGInterpolation(SolverOptions.CLASSIC_INTERPOLATION_WITH_FF_COUPLING)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_PoissonSqueezedX(self):
        x=self.domain.getX().copy()
        x[0]*=0.5
        self.domain.setX(x)
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(1,Solution(self.domain))
        g_ex=Vector(0.,Solution(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)


   def test_Poisson2(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_Poisson3(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_Poisson4(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)
   def test_Poisson2Classic(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        pde.getSolverOptions().setAMGInterpolation(SolverOptions.CLASSIC_INTERPOLATION)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_Poisson3Classic(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        pde.getSolverOptions().setAMGInterpolation(SolverOptions.CLASSIC_INTERPOLATION)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_Poisson4Classic(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        pde.getSolverOptions().setAMGInterpolation(SolverOptions.CLASSIC_INTERPOLATION)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_WeakCoupled2(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1

        Y=Data(0.,(2,),Function(self.domain))
        a=-1./2.*0.01
        Y[0]=    u_ex[0] + a * u_ex[1]
        Y[1]=a * u_ex[0] +     u_ex[1]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=2)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=kronecker(2)+a*(numpy.ones((2,2))-kronecker(2)),
                     Y=Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_WeakCoupled3(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1

        a=-1./3.*0.01
        Y=Data(0.,(3,),Function(self.domain))
        Y[0]=     u_ex[0]+ a * u_ex[1]+ a * u_ex[2]
        Y[1]= a * u_ex[0]+     u_ex[1]+ a * u_ex[2]
        Y[2]= a * u_ex[0]+ a * u_ex[1]+     u_ex[2]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=3)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=kronecker(3)+a*(numpy.ones((3,3))-kronecker(3)),
                     Y=Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_WeakCoupled4(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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

        a=-1./4.*0.01

        Y=Data(0.,(4,),Function(self.domain))
        Y[0]=     u_ex[0]+ a * u_ex[1]+ a * u_ex[2]+ a * u_ex[3]
        Y[1]= a * u_ex[0]+     u_ex[1]+ a * u_ex[2]+ a * u_ex[3]
        Y[2]= a * u_ex[0]+ a * u_ex[1]+     u_ex[2]+ a * u_ex[3]
        Y[3]= a * u_ex[0]+ a * u_ex[1]+ a * u_ex[2]+     u_ex[3]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=4)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=kronecker(4)+a*(numpy.ones((4,4))-kronecker(4)),
                     Y=Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_StrongCoupled2(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(2,),Solution(self.domain))
        g_ex=Data(0.,(2,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(2,self.domain.getDim(),2,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1

        Y=Data(0.,(2,),Function(self.domain))
        a=-1./2.*0.99
        Y[0]=    u_ex[0] + a * u_ex[1]
        Y[1]=a * u_ex[0] +     u_ex[1]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=2)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=kronecker(2)+a*(numpy.ones((2,2))-kronecker(2)),
                     Y=Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)

        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_StrongCoupled3(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(3,),Solution(self.domain))
        g_ex=Data(0.,(3,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(3,self.domain.getDim(),3,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
           u_ex[0]+= 1*(i+1)*x[i]
           g_ex[0,i]=1*(i+1)
           u_ex[1]+= 2*(i+1)*x[i]
           g_ex[1,i]=2*(i+1)
           u_ex[2]+= 3*(i+1)*x[i]
           g_ex[2,i]=3*(i+1)
           A[0,i,0,i]=1
           A[1,i,1,i]=1
           A[2,i,2,i]=1

        a=-1./3.*0.99
        Y=Data(0.,(3,),Function(self.domain))
        Y[0]=     u_ex[0]+ a * u_ex[1]+ a * u_ex[2]
        Y[1]= a * u_ex[0]+     u_ex[1]+ a * u_ex[2]
        Y[2]= a * u_ex[0]+ a * u_ex[1]+     u_ex[2]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=3)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=kronecker(3)+a*(numpy.ones((3,3))-kronecker(3)),
                     Y=Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)

   def test_StrongCoupled4(self):
        x=Solution(self.domain).getX()
        # --- set exact solution ----
        u_ex=Data(1.,(4,),Solution(self.domain))
        g_ex=Data(0.,(4,self.domain.getDim()), Solution(self.domain))
        A=Data(0.,(4,self.domain.getDim(),4,self.domain.getDim()), Function(self.domain))
        for i in range(self.domain.getDim()):
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

        a=-1./4.*0.99

        Y=Data(0.,(4,),Function(self.domain))
        Y[0]=     u_ex[0]+ a * u_ex[1]+ a * u_ex[2]+ a * u_ex[3]
        Y[1]= a * u_ex[0]+     u_ex[1]+ a * u_ex[2]+ a * u_ex[3]
        Y[2]= a * u_ex[0]+ a * u_ex[1]+     u_ex[2]+ a * u_ex[3]
        Y[3]= a * u_ex[0]+ a * u_ex[1]+ a * u_ex[2]+     u_ex[3]
        # create PDE:
        pde=LinearPDE(self.domain,numEquations=4)
        pde.setSymmetryOn()
        mask=whereZero(x[0])*[1,1,1,1]
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=A,
                     D=kronecker(4)+a*(numpy.ones((4,4))-kronecker(4)),
                     Y=Y,
                     y=matrixmult(g_ex,self.domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)


   def test_Square(self):
        # PDE constants
        h1 = 0.5
        h2 = 0.5
        beta = 1.
        alpha1 = 2.
        alpha2 = 1.

        # domain masks and domain-specific constants
        x = Solution(self.domain).getX(); x0 = x[0]; x1 = x[1]
        omega2 = wherePositive(x0-h1)*wherePositive(x1-h2)
        omega1 = 1-omega2
        ratio = alpha1/alpha2
        alpha = alpha1*omega1 + alpha2*omega2

        # --- set exact solution ----
        a1 = 1.
        d1 = 1.
        b1 = -d1*h2
        c1 = -d1*h1

        a2 = a1 - (1.-ratio)*d1*h1*h2
        b2 = ratio*b1
        c2 = ratio*c1
        d2 = ratio*d1

        u_ex = omega1*(a1 + b1*x0 + c1*x1 + d1*x0*x1) + \
                omega2*(a2 + b2*x0 + c2*x1 + d2*x0*x1)

        # create PDE:
        pde = LinearPDE(self.domain,numEquations=1)

        # set the value to that of the solution on the boundary
        q = whereZero(x0) + whereZero(x1) + \
            whereZero(sup(x0)-x0) + whereZero(sup(x1)-x1)
        pde.setValue(q=q,r=u_ex)

        # create X points in the centre of the grid elements
        xe = Function(self.domain).getX()
        x0 = xe[0]
        x1 = xe[1]

        # redefine omega so that apha is more precise on the diagonal (?)
        omega2 = wherePositive(x0-h1)*wherePositive(x1-h2)
        omega1 = 1-omega2
        ratio = alpha1/alpha2
        alpha = alpha1*omega1 + alpha2*omega2

        # set up PDE coefficients
        pde.setValue(A=alpha*kronecker(self.domain), D=beta, Y=beta*u_ex)
        pde.setSymmetryOn()

        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(self.SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        if MIN_MATRIX_SIZE!= None: pde.getSolverOptions().setMinCoarseMatrixSize(MIN_MATRIX_SIZE)
        if MIN_SPARSITY!=None: pde.getSolverOptions().setMinCoarseMatrixSparsity(MIN_SPARSITY)
        if MAX_LEVEL!=None: pde.getSolverOptions().setLevelMax(MAX_LEVEL)
        u = pde.getSolution()

        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.assertTrue(error<self.RES_TOL, "solution error %s is too big."%error)


class Test_AMGOnFinleyHex2DOrder1(AMG):
   RES_TOL=5.e-7
   SOLVER_TOL=1.e-8
   def setUp(self):
        NE=int(float(NE_TOTAL)**(1./2.)+0.5)
        self.domain = Rectangle(NE,NE,1, optimize=OPTIMIZE)
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

