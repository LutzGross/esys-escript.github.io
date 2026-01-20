
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

NUM_THREADS=8
import os
TEST_STR="timing: per iteration step:"
REPEAT=10
HEADER="""from esys.escript import *
from esys.finley import Rectangle,Brick 
from esys.escript.linearPDEs import LinearPDE 
SOLVER_TOL=1.e-2
REL_TOL=1.
OPTIMIZE=False 
SOLVER_VERBOSE=True
FAC_DIAG=1.
FAC_OFFDIAG=-0.4

setNumberOfThreads(%d)
"""

DOM_2_1="dom=Rectangle(NE,NE,order=1, useFullElementOrder=False,optimize=OPTIMIZE)"
DOM_2_2="dom=Rectangle(NE,NE,order=2, useFullElementOrder=False,optimize=OPTIMIZE)"
DOM_3_1="dom=Brick(NE,NE,NE,order=1, useFullElementOrder=True,optimize=OPTIMIZE)"
DOM_3_2="dom=Brick(NE,NE,NE,order=2, useFullElementOrder=True,optimize=OPTIMIZE)"

TEST_2_s="""x=Solution(dom).getX()
u_ex=Scalar(0,Solution(dom))
u_ex=1.+2.*x[0]+3.*x[1]
g_ex=Data(0.,(2,),Solution(dom))
g_ex[0]=2.
g_ex[1]=3.
pde=LinearPDE(dom,numEquations=1)
mask=whereZero(x[0])
pde.setValue(r=u_ex,q=mask)
pde.setValue(A=kronecker(2),y=inner(g_ex,dom.getNormal()))
"""
TEST_2_v="""x=Solution(dom).getX()
x=Solution(dom).getX()
u_ex=Vector(0,Solution(dom))
u_ex[0]=1.+2.*x[0]+3.*x[1]
u_ex[1]=-1.+3.*x[0]+2.*x[1]
g_ex=Data(0.,(2,2),Solution(dom))
g_ex[0,0]=2.
g_ex[0,1]=3.
g_ex[1,0]=3.
g_ex[1,1]=2.
pde=LinearPDE(dom,numEquations=2)
mask=whereZero(x[0])
pde.setValue(r=u_ex,q=mask*numarray.ones(2,))
A=Tensor4(0,Function(dom))
A[0,:,0,:]=kronecker(2)
A[1,:,1,:]=kronecker(2)
Y=Vector(0.,Function(dom))
Y[0]=u_ex[0]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG
Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG
pde.setValue(A=A, D=kronecker(2)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((2,2))*FAC_OFFDIAG, Y=Y, y=matrixmult(g_ex,dom.getNormal()))
"""

TEST_3_s="""x=Solution(dom).getX()
u_ex=1.+2.*x[0]+3.*x[1]+4.*x[2]
g_ex=Data(0.,(3,),Solution(dom))
g_ex[0]=2.
g_ex[1]=3.
g_ex[2]=4.
pde=LinearPDE(dom,numEquations=1)
mask=whereZero(x[0])
pde.setValue(r=u_ex,q=mask)
pde.setValue(A=kronecker(3),y=inner(g_ex,dom.getNormal()))
"""

TEST_3_v="""x=Solution(dom).getX()
u_ex=Vector(0,Solution(dom))
u_ex[0]=1.+2.*x[0]+3.*x[1]+4.*x[2]
u_ex[1]=-1.+4.*x[0]+1.*x[1]-2.*x[2]
u_ex[2]=5.+8.*x[0]+4.*x[1]+5.*x[2]
g_ex=Data(0.,(3,3),Solution(dom))
g_ex[0,0]=2.
g_ex[0,1]=3.
g_ex[0,2]=4.
g_ex[1,0]=4.
g_ex[1,1]=1.
g_ex[1,2]=-2.
g_ex[2,0]=8.
g_ex[2,1]=4.
g_ex[2,2]=5.
pde=LinearPDE(dom,numEquations=3)
mask=whereZero(x[0])
pde.setValue(r=u_ex,q=mask*numarray.ones(3,))
A=Tensor4(0,Function(dom))
A[0,:,0,:]=kronecker(3)
A[1,:,1,:]=kronecker(3)
A[2,:,2,:]=kronecker(3)
Y=Vector(0.,Function(dom))
Y[0]=u_ex[0]*FAC_DIAG+u_ex[2]*FAC_OFFDIAG+u_ex[1]*FAC_OFFDIAG
Y[1]=u_ex[1]*FAC_DIAG+u_ex[0]*FAC_OFFDIAG+u_ex[2]*FAC_OFFDIAG
Y[2]=u_ex[2]*FAC_DIAG+u_ex[1]*FAC_OFFDIAG+u_ex[0]*FAC_OFFDIAG
pde.setValue(A=A,
D=kronecker(3)*(FAC_DIAG-FAC_OFFDIAG)+numarray.ones((3,3))*FAC_OFFDIAG,
Y=Y,
y=matrixmult(g_ex,dom.getNormal()))
"""

SOLVE_AND_TEST="""pde.setTolerance(SOLVER_TOL)
pde.setSolverMethod(pde.PCG,pde.JACOBI)
pde.setSolverPackage(pde.PASO)
u=pde.getSolution(verbose=SOLVER_VERBOSE)
error=Lsup(u-u_ex)/Lsup(u_ex)
if error>REL_TOL*Lsup(u_ex): raise RuntimeError("solution error %s is too big."%error)
"""


#for n in [10000, 50000, 100000]:
for n in [100000]:
# for n in [1000, 10000]:
 #for prop in [ (1,2), (2,2), (1,3), (2,3) ]:
 for prop in [ (1,2), (1,3) ]:
   for tp in [ "s", "v" ]:
      # create code:
      prog=HEADER%NUM_THREADS
      dim=prop[1]
      if isinstance(prop[0], int):
          o=prop[0]
          if tp=="s": 
                q=1
          else:
                q=dim
          NE=int(float(n/q-1)**(1./dim)/o+0.5)
          prog+="NE=%d\n"%NE
          if dim==2:
              if o==1:
                 prog+=DOM_2_1
              else:
                 prog+=DOM_2_2
          else:
              if o==1:
                 prog+=DOM_3_1
              else:
                 prog+=DOM_3_2
          prog+="\n"
      if dim==2:
        if tp =="s":
           prog+=TEST_2_s
        else:
           prog+=TEST_2_v
      else:
        if tp =="s":
           prog+=TEST_3_s
        else:
           prog+=TEST_3_v
      print("l= %d, dim= %d, type=%s, order=%s"%(q*(o*NE+1)**dim,dim,tp,o))
    
      prog+=SOLVE_AND_TEST 
      # run code:
      print(prog, file=file("__prog","w"))
      # activate for dynamic
      # for CHUNK in [1,10,100,1000,10000, 100000]:
      #   for CHUNK_PCG in [1,10,100,1000,10000, 100000]:
      # activate for static
      for CHUNK in [-1]:
       for CHUNK_PCG in [-1]:
        if CHUNK*NUM_THREADS <= n and CHUNK_PCG*NUM_THREADS <=n:
         time_per_iter=0
         for i in range(REPEAT):
            os.system("export OMP_NUM_THREADS=%d;export PASO_CHUNK_SIZE_MVM=%d; export PASO_CHUNK_SIZE_PCG=%d; python __prog > __out;"%(NUM_THREADS,CHUNK,CHUNK_PCG))
            out=file("__out","r").read()
            for i in out.split("\n"):
               if i.startswith(TEST_STR): time_per_iter+=float(i[len(TEST_STR):-3].strip())
         print(CHUNK,CHUNK_PCG,time_per_iter/REPEAT)

