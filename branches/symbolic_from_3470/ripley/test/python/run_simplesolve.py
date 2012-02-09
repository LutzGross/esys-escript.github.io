
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
Test suite for the linearPDE  and pdetools test on ripley

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
from esys.ripley import Rectangle,Brick
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import numpy
SOLVER_VERBOSE=False 

try:
     RIPLEY_TEST_DATA=os.environ['RIPLEY_TEST_DATA']
except KeyError:
     RIPLEY_TEST_DATA='.'

SOLVER_TOL=1.e-8
REL_TOL=1.e-6

FAC_DIAG=1.
FAC_OFFDIAG=-0.4

# number of elements in the spatial directions
NE0=12
NE1=12
NE2=8
mpiSize=getMPISizeWorld()
for x in [int(sqrt(mpiSize)),2,3,5,7,1]:
    NX=x
    NY=mpiSize/x
    if NX*NY == mpiSize:
        break

for x in [(int(mpiSize**(1/3.)),int(mpiSize**(1/3.))),(2,3),(2,2),(1,2),(1,1)]:
    NXb=x[0]
    NYb=x[1]
    NZb=mpiSize/(x[0]*x[1])
    if NXb*NYb*NZb == mpiSize:
        break

class SimpleSolve_Rectangle_SinglePDE_Paso_BICGSTAB_Jacobi(unittest.TestCase):
     def test_solve(self):
        # Tell about how many MPI CPUs and OpenMP threads
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(0,Solution(domain))
        u_ex=1.+2.*x[0]+3.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        # -------- test gradient --------------------------------
        g=grad(u_ex)
        self.assertTrue(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.BICGSTAB)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Rectangle_SinglePDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(0,Solution(domain))
        u_ex=1.+2.*x[0]+3.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        # -------- test gradient --------------------------------
        g=grad(u_ex)
        self.assertTrue(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Rectangle_SystemPDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]
        u_ex[1]=-1.+3.*x[0]+2.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,2),Solution(domain))
        g_ex[0,0]=2.
        g_ex[0,1]=3.
        g_ex[1,0]=3.
        g_ex[1,1]=2.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numpy.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((2,2))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_SinglePDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        x=Solution(domain).getX()
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        g_ex[2]=4.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_SystemPDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]
        u_ex[1]=-1.+4.*x[0]+1.*x[1]-2.*x[2]
        u_ex[2]=5.+8.*x[0]+4.*x[1]+5.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,3),Solution(domain))
        g_ex[0,0]=2.
        g_ex[0,1]=3.
        g_ex[0,2]=4.
        g_ex[1,0]=4.
        g_ex[1,1]=1.
        g_ex[1,2]=-2.
        g_ex[2,0]=8.
        g_ex[2,1]=4.
        g_ex[2,2]=5.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numpy.ones(3,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((3,3))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_SinglePDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(0,Solution(domain))
        u_ex=1.+2.*x[0]+3.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        # -------- test gradient --------------------------------
        g=grad(u_ex)
        self.assertTrue(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_SystemPDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]
        u_ex[1]=-1.+3.*x[0]+2.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,2),Solution(domain))
        g_ex[0,0]=2.
        g_ex[0,1]=3.
        g_ex[1,0]=3.
        g_ex[1,1]=2.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numpy.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((2,2))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_SinglePDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        x=Solution(domain).getX()
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        g_ex[2]=4.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_SystemPDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]
        u_ex[1]=-1.+4.*x[0]+1.*x[1]-2.*x[2]
        u_ex[2]=5.+8.*x[0]+4.*x[1]+5.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,3),Solution(domain))
        g_ex[0,0]=2.
        g_ex[0,1]=3.
        g_ex[0,2]=4.
        g_ex[1,0]=4.
        g_ex[1,1]=1.
        g_ex[1,2]=-2.
        g_ex[2,0]=8.
        g_ex[2,1]=4.
        g_ex[2,2]=5.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numpy.ones(3,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((3,3))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.TFQMR)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Rectangle_SinglePDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Scalar(0,Solution(domain))
        u_ex=1.+2.*x[0]+3.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        # -------- test gradient --------------------------------
        g=grad(u_ex)
        self.assertTrue(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_SystemPDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(n0=NE0*NX-1, n1=NE1*NY-1, d0=NX, d1=NY)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]
        u_ex[1]=-1.+3.*x[0]+2.*x[1]
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,2),Solution(domain))
        g_ex[0,0]=2.
        g_ex[0,1]=3.
        g_ex[1,0]=3.
        g_ex[1,1]=2.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numpy.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((2,2))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_SinglePDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        x=Solution(domain).getX()
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        g_ex[2]=4.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_SystemPDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(n0=NE0*NXb-1, n1=NE1*NYb-1, n2=NE2*NZb-1, d0=NXb, d1=NYb, d2=NZb)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]
        u_ex[1]=-1.+4.*x[0]+1.*x[1]-2.*x[2]
        u_ex[2]=5.+8.*x[0]+4.*x[1]+5.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,3),Solution(domain))
        g_ex[0,0]=2.
        g_ex[0,1]=3.
        g_ex[0,2]=4.
        g_ex[1,0]=4.
        g_ex[1,1]=1.
        g_ex[1,2]=-2.
        g_ex[2,0]=8.
        g_ex[2,1]=4.
        g_ex[2,2]=5.
        # -------- test gradient --------------------------------
        self.assertTrue(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numpy.ones(3,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numpy.ones((3,3))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.getSolverOptions().setTolerance(SOLVER_TOL)
        pde.getSolverOptions().setSolverMethod(SolverOptions.MINRES)
        pde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
        pde.getSolverOptions().setPackage(SolverOptions.PASO)
        pde.getSolverOptions().setVerbosity(SOLVER_VERBOSE)
        u=pde.getSolution()
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)
        self.assertTrue(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SinglePDE_Paso_BICGSTAB_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SinglePDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SystemPDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_SinglePDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_SystemPDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SinglePDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SystemPDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_SinglePDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_SystemPDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_SystemPDE_Paso_MINRES_Jacobi))
 
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)
