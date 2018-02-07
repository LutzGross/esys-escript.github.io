
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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
#
# this script is testing block solvers for PDEs
#
#
#    - u_{j,ii} + b*u_j+ a*sum_{k<>j}  (u_j-u_k) = F_j
#
#  where a controls the degree of coupling and b the degree of diagonal dominance.
#  a and b may have any value.
#  
#  The domain needs to be a unit square or cube with any type of mesh
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
import esys.dudley as dudley

TOL=1.e-8

def runTest(dom, n=1, a=0, b=0):
   print("================== TEST : n= %s a=%s b=%s ================"%(n,a,b))
   DIM=dom.getDim()
   normal=dom.getNormal()
   mypde=LinearPDE(dom, numEquations=n, numSolutions=n)
   x=dom.getX()
   A=mypde.createCoefficient("A")
   D=mypde.createCoefficient("D")
   Y=mypde.createCoefficient("Y")
   y=mypde.createCoefficient("y")
   q=mypde.createCoefficient("q")
   if n==1:
      u_ref=Scalar(0.,Solution(dom))
      for j in range(DIM):
         q+=whereZero(x[j])
         A[j,j]=1
      y+=sin(sqrt(2.))*normal[0]
      Y+=b*x[0]*sin(sqrt(2.))
      D+=b
      u_ref=x[0]*sin(sqrt(2.))
   else:
      u_ref=Data(0.,(n,),Solution(dom))
      for i in range(n):
          for j in range(DIM):
             q[i]+=whereZero(x[j])
             A[i,j,i,j]=1
             if j == i%DIM: y[i]+=sin(i+sqrt(2.))*normal[j]
          for j in range(n):
                if i==j:
                   D[i,i]=b+(n-1)*a
                   Y[i]+=b*x[i%DIM]*sin(i+sqrt(2.))
                else:
                   D[i,j]=-a
                   Y[i]+=a*(x[i%DIM]*sin(i+sqrt(2.))-x[j%DIM]*sin(j+sqrt(2.)))
          u_ref[i]=x[i%DIM]*sin(i+sqrt(2.))
          
#    - u_{j,ii} + b*u_i+ a*sum_{k<>j}  (u_i-u_k) = F_j
          
   mypde.setValue(A=A, D=D, y=y, q=q, r=u_ref, Y=Y)
   mypde.getSolverOptions().setVerbosityOn()
   mypde.getSolverOptions().setTolerance(TOL)
   mypde.setSymmetryOn()
   u=mypde.getSolution()
   error=Lsup(u-u_ref)/Lsup(u_ref)
   print("error = ",error)
   if error > 10*TOL: print("XXXXXXXXXX Error to large ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

domain=dudley.Rectangle(10,20,1,l0=1.,l1=1.0)
# or Brick or ReadMesh
runTest(dom=domain, n=1, b=0)
runTest(dom=domain, n=1, b=5)
runTest(dom=domain, n=1, b=50)

runTest(dom=domain, n=2, a=0, b=0)
runTest(dom=domain, n=2, a=5, b=0)
runTest(dom=domain, n=2, a=50, b=0)

runTest(dom=domain, n=2, a=0, b=10)
runTest(dom=domain, n=2, a=5, b=10)
runTest(dom=domain, n=2, a=50, b=10)

runTest(dom=domain, n=3, a=0, b=0)
runTest(dom=domain, n=3, a=5, b=0)
runTest(dom=domain, n=3, a=50, b=0)

runTest(dom=domain, n=3, a=0, b=10)
runTest(dom=domain, n=3, a=5, b=10)
runTest(dom=domain, n=3, a=50, b=10)
