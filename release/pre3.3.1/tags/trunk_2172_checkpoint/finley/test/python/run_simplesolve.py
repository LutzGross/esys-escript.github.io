
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
Test suite for the linearPDE  and pdetools test on finley

@remark:

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import unittest, sys

from esys.escript import *
from esys.finley import Rectangle,Brick
from esys.escript.linearPDEs import LinearPDE
OPTIMIZE=True
SOLVER_VERBOSE=False 
# setNumberOfThreads(2)

try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=FINLEY_TEST_DATA+"/data_meshes/"

# number of elements in the spatial directions
NE0=8
NE1=10
NE2=12

NE0=12
NE1=12
NE2=8

SOLVER_TOL=1.e-8
REL_TOL=1.e-6

FAC_DIAG=1.
FAC_OFFDIAG=-0.4


class SimpleSolve_Rectangle_Order1_SinglePDE_Paso_BICGSTAB_Jacobi(unittest.TestCase):
     def test_solve(self):
	# Tell about how many MPI CPUs and OpenMP threads
        domain=Rectangle(NE0,NE1,1, optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.BICGSTAB,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Rectangle_Order1_SinglePDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,1, optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Rectangle_Order1_SystemPDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,1,optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Rectangle_Order2_SinglePDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,2,l0=1.,l1=1,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.+8.*x[0]+5.*x[1]
        g_ex[1]=3.+5.*x[0]+12.*x[1]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()),Y=-20.)
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Rectangle_Order2_SystemPDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
        u_ex[1]=-1.+4.*x[0]+2.*x[1]+1.*x[0]**2+6.*x[1]*x[0]+4.*x[1]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,2),Solution(domain))
        g_ex[0,0]=2.+8.*x[0]+5.*x[1]
        g_ex[0,1]=3.+5.*x[0]+12.*x[1]
        g_ex[1,0]=4.+2.*x[0]+6.*x[1]
        g_ex[1,1]=2.+6.*x[0]+8.*x[1]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG,
                     Y=Y-[20.,10.],
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_Order1_SinglePDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,1,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        g_ex[2]=4.
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_Order1_SystemPDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,1,optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_Order2_SinglePDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]+6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.+6.*x[1]+8.*x[2]+18.*x[0]
        g_ex[1]=3.+6.*x[0]+7.*x[2]+20.*x[1]
        g_ex[2]=4.+7.*x[1]+8.*x[0]+22.*x[2]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()),Y=-60.)
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
class SimpleSolve_Brick_Order2_SystemPDE_Paso_PCG_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]+6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        u_ex[1]=2.+4.*x[0]+1.*x[1]-6.*x[2]+3.*x[0]*x[1]+2.*x[1]*x[2]-8.*x[2]*x[0]-2.*x[0]**2+7.*x[1]**2+5.*x[2]**2
        u_ex[2]=-2.+7.*x[0]+9.*x[1]+2*x[2]-6.*x[0]*x[1]+8.*x[1]*x[2]+2.*x[2]*x[0]+2.*x[0]**2+8.*x[1]**2+1.*x[2]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,3),Solution(domain))
        g_ex[0,0]=2.+6.*x[1]+8.*x[2]+18.*x[0]
        g_ex[0,1]=3.+6.*x[0]+7.*x[2]+20.*x[1]
        g_ex[0,2]=4.+7.*x[1]+8.*x[0]+22.*x[2]
        g_ex[1,0]=4.+3.*x[1]-8.*x[2]-4.*x[0]
        g_ex[1,1]=1+3.*x[0]+2.*x[2]+14.*x[1]
        g_ex[1,2]=-6.+2.*x[1]-8.*x[0]+10.*x[2]
        g_ex[2,0]=7.-6.*x[1]+2.*x[2]+4.*x[0]
        g_ex[2,1]=9.-6.*x[0]+8.*x[2]+16.*x[1]
        g_ex[2,2]=2+8.*x[1]+2.*x[0]+2.*x[2]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
                     Y=Y-numarray.array([60.,20.,22.]),
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.PCG,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_Order1_SinglePDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,1, optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_Order2_SinglePDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,2,l0=1.,l1=1,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.+8.*x[0]+5.*x[1]
        g_ex[1]=3.+5.*x[0]+12.*x[1]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()),Y=-20.)
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Rectangle_Order1_SystemPDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,1,optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_Order2_SystemPDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
        u_ex[1]=-1.+4.*x[0]+2.*x[1]+1.*x[0]**2+6.*x[1]*x[0]+4.*x[1]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,2),Solution(domain))
        g_ex[0,0]=2.+8.*x[0]+5.*x[1]
        g_ex[0,1]=3.+5.*x[0]+12.*x[1]
        g_ex[1,0]=4.+2.*x[0]+6.*x[1]
        g_ex[1,1]=2.+6.*x[0]+8.*x[1]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG,
                     Y=Y-[20.,10.],
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Brick_Order1_SinglePDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,1,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        g_ex[2]=4.
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_Order1_SystemPDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,1,optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_Order2_SinglePDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]+6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.+6.*x[1]+8.*x[2]+18.*x[0]
        g_ex[1]=3.+6.*x[0]+7.*x[2]+20.*x[1]
        g_ex[2]=4.+7.*x[1]+8.*x[0]+22.*x[2]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()),Y=-60.)
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_Order2_SystemPDE_Paso_TFQMR_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]+6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        u_ex[1]=2.+4.*x[0]+1.*x[1]-6.*x[2]+3.*x[0]*x[1]+2.*x[1]*x[2]-8.*x[2]*x[0]-2.*x[0]**2+7.*x[1]**2+5.*x[2]**2
        u_ex[2]=-2.+7.*x[0]+9.*x[1]+2*x[2]-6.*x[0]*x[1]+8.*x[1]*x[2]+2.*x[2]*x[0]+2.*x[0]**2+8.*x[1]**2+1.*x[2]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,3),Solution(domain))
        g_ex[0,0]=2.+6.*x[1]+8.*x[2]+18.*x[0]
        g_ex[0,1]=3.+6.*x[0]+7.*x[2]+20.*x[1]
        g_ex[0,2]=4.+7.*x[1]+8.*x[0]+22.*x[2]
        g_ex[1,0]=4.+3.*x[1]-8.*x[2]-4.*x[0]
        g_ex[1,1]=1+3.*x[0]+2.*x[2]+14.*x[1]
        g_ex[1,2]=-6.+2.*x[1]-8.*x[0]+10.*x[2]
        g_ex[2,0]=7.-6.*x[1]+2.*x[2]+4.*x[0]
        g_ex[2,1]=9.-6.*x[0]+8.*x[2]+16.*x[1]
        g_ex[2,2]=2+8.*x[1]+2.*x[0]+2.*x[2]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
                     Y=Y-numarray.array([60.,20.,22.]),
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.TFQMR,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)


class SimpleSolve_Rectangle_Order1_SinglePDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,1, optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-g)<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Rectangle_Order2_SinglePDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,2,l0=1.,l1=1,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,),Solution(domain))
        g_ex[0]=2.+8.*x[0]+5.*x[1]
        g_ex[1]=3.+5.*x[0]+12.*x[1]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(2),y=inner(g_ex,domain.getNormal()),Y=-20.)
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Rectangle_Order1_SystemPDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,1,optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Rectangle_Order2_SystemPDE_Paso_MINRES_RILU(unittest.TestCase):
     def test_solve(self):
        domain=Rectangle(NE0,NE1,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[0]**2+5.*x[1]*x[0]+6.*x[1]**2
        u_ex[1]=-1.+4.*x[0]+2.*x[1]+1.*x[0]**2+6.*x[1]*x[0]+4.*x[1]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(2,2),Solution(domain))
        g_ex[0,0]=2.+8.*x[0]+5.*x[1]
        g_ex[0,1]=3.+5.*x[0]+12.*x[1]
        g_ex[1,0]=4.+2.*x[0]+6.*x[1]
        g_ex[1,1]=2.+6.*x[0]+8.*x[1]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=2)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(2)
        A[1,:,1,:]=kronecker(2)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG,
                     Y=Y-[20.,10.],
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.RILU)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

class SimpleSolve_Brick_Order1_SinglePDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,1,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.
        g_ex[1]=3.
        g_ex[2]=4.
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_Order1_SystemPDE_Paso_MINRES_Jacobi(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,1,optimize=OPTIMIZE)
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
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
                     Y=Y,
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.JACOBI)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_Order2_SinglePDE_Paso_MINRES_RILU(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]+6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,),Solution(domain))
        g_ex[0]=2.+6.*x[1]+8.*x[2]+18.*x[0]
        g_ex[1]=3.+6.*x[0]+7.*x[2]+20.*x[1]
        g_ex[2]=4.+7.*x[1]+8.*x[0]+22.*x[2]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=1)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask)
        pde.setValue(A=kronecker(3),y=inner(g_ex,domain.getNormal()),Y=-60.)
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.RILU)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)
        
class SimpleSolve_Brick_Order2_SystemPDE_Paso_MINRES_RILU(unittest.TestCase):
     def test_solve(self):
        domain=Brick(NE0,NE1,NE2,2,optimize=OPTIMIZE)
        x=Solution(domain).getX()
        # --- set exact solution ----
        u_ex=Vector(0,Solution(domain))
        u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]+6.*x[0]*x[1]+7.*x[1]*x[2]+8.*x[2]*x[0]+9.*x[0]**2+10.*x[1]**2+11.*x[2]**2
        u_ex[1]=2.+4.*x[0]+1.*x[1]-6.*x[2]+3.*x[0]*x[1]+2.*x[1]*x[2]-8.*x[2]*x[0]-2.*x[0]**2+7.*x[1]**2+5.*x[2]**2
        u_ex[2]=-2.+7.*x[0]+9.*x[1]+2*x[2]-6.*x[0]*x[1]+8.*x[1]*x[2]+2.*x[2]*x[0]+2.*x[0]**2+8.*x[1]**2+1.*x[2]**2
        # --- set exact gradient -----------
        g_ex=Data(0.,(3,3),Solution(domain))
        g_ex[0,0]=2.+6.*x[1]+8.*x[2]+18.*x[0]
        g_ex[0,1]=3.+6.*x[0]+7.*x[2]+20.*x[1]
        g_ex[0,2]=4.+7.*x[1]+8.*x[0]+22.*x[2]
        g_ex[1,0]=4.+3.*x[1]-8.*x[2]-4.*x[0]
        g_ex[1,1]=1+3.*x[0]+2.*x[2]+14.*x[1]
        g_ex[1,2]=-6.+2.*x[1]-8.*x[0]+10.*x[2]
        g_ex[2,0]=7.-6.*x[1]+2.*x[2]+4.*x[0]
        g_ex[2,1]=9.-6.*x[0]+8.*x[2]+16.*x[1]
        g_ex[2,2]=2+8.*x[1]+2.*x[0]+2.*x[2]
        # -------- test gradient --------------------------------
        self.failUnless(Lsup(g_ex-grad(u_ex))<REL_TOL*Lsup(g_ex))
        # -------- set-up PDE ----------------------------------- 
        pde=LinearPDE(domain,numEquations=3)
        mask=whereZero(x[0])
        pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
        Y=Vector(0.,Function(domain))
        Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
        Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
        Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
        A=Tensor4(0,Function(domain))
        A[0,:,0,:]=kronecker(3)
        A[1,:,1,:]=kronecker(3)
        A[2,:,2,:]=kronecker(3)
        pde.setValue(A=A,
                     D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
                     Y=Y-numarray.array([60.,20.,22.]),
                     y=matrixmult(g_ex,domain.getNormal()))
        # -------- get the solution ---------------------------
        pde.setTolerance(SOLVER_TOL)
        pde.setSolverMethod(pde.MINRES,pde.RILU)
        pde.setSolverPackage(pde.PASO)
        u=pde.getSolution(verbose=SOLVER_VERBOSE)
        # -------- test the solution ---------------------------
        error=Lsup(u-u_ex)/Lsup(u_ex)
        self.failUnless(error<REL_TOL*Lsup(u_ex), "solution error %s is too big."%error)

if __name__ == '__main__':
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SinglePDE_Paso_BICGSTAB_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SinglePDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SystemPDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order2_SinglePDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order2_SystemPDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order1_SinglePDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order1_SystemPDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order2_SinglePDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order2_SystemPDE_Paso_PCG_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SinglePDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order2_SinglePDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SystemPDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order2_SystemPDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order1_SinglePDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order1_SystemPDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order2_SinglePDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order2_SystemPDE_Paso_TFQMR_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order2_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order1_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Rectangle_Order2_SystemPDE_Paso_MINRES_RILU))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order1_SinglePDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order1_SystemPDE_Paso_MINRES_Jacobi))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order2_SinglePDE_Paso_MINRES_RILU))
   suite.addTest(unittest.makeSuite(SimpleSolve_Brick_Order2_SystemPDE_Paso_MINRES_RILU))
 
   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)
