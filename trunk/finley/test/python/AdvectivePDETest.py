
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

# Test for the AdvectivePDE class
#
#  for a single equation the test problem is
#
#   -(K_{ij}u_{,j})_{,i} - (w_i u)_{,i} + v_j u_{,j} =0 
#
#   + constraints on the surface
#
#  for system of two equation the test problem is
#
#   -(K_{milj}u_{l,j})_{,i} - (w_{mil} u_l)_{,i} + v_{mlj} u_{l,j} =0 
#
#   + constraints on the surface
#
#   K,w and v are constant (we will set v=0 or w=0)
#
#   the test solution is  u(x)=e^{z_i*x_i}  and u_l(x)=e^{z_{li}*x_i}
#
#   an easy caculation shows that    
#    
#     z_i*K_{ij}*z_j=(v_i-w_i)*z_i and  z_{li}*K_{milj}*z_{lj}=(v_{mjl}-w_{mlj})*z_{lj}
#
#   obviously one can choose: v_i-w_i=K_{ji}z_j and v_{mjl}-w_{mlj}=z_{li}*K_{milj} (no summation over l)
#

from esys.escript import *
from esys.escript.linearPDEs import AdvectivePDE,LinearPDE
from esys import finley
from random import random

def printError(u,u_ex):
    if u.getRank()==0:
       out=" error = %e range = [%e:%e] [%e:%e]"%(Lsup(u-u_ex)/Lsup(u_ex),sup(u),inf(u),sup(u_ex),inf(u_ex))
    else:
       out="\n"
       for i in range(u.getShape()[0]):
          out+="     %d error = %e range = [%e:%e] [%e:%e]\n"%(i,Lsup(u[i]-u_ex[i])/Lsup(u_ex[i]),sup(u[i]),inf(u[i]),sup(u_ex[i]),inf(u_ex[i]))
    return out
    
    
def makeRandomFloats(n,val_low=0.,val_up=1.):
    out=[]
    for i in range(n):
         out.append((val_up-val_low)*random()+val_low)
    return out

def makeRandomFloatMatrix(m,n,val_low=0.,val_up=1.):
    out=[]
    for i in range(m):
         out.append(makeRandomFloats(n,val_low,val_up))
    return out

def makeRandomFloatTensor(l,k,m,n,val_low=0.,val_up=1.):
    out=[]
    for j in range(l):
         out2=[]
         for i in range(k): out2.append(makeRandomFloatMatrix(m,n,val_low,val_up))
         out.append(out2)
    return out

ne=20
# for d in [2,3]:
for d in [3]:
   # create domain:
   if d==2:
     mydomain=finley.Rectangle(ne,ne,1)
     x=mydomain.getX()
     msk=whereZero(x[0])+whereZero(x[0]-1.)+whereZero(x[1])+whereZero(x[1]-1.)
   else:
     mydomain=finley.Brick(ne,ne,ne,1)
     x=mydomain.getX()
     msk=whereZero(x[0])+whereZero(x[0]-1.)+whereZero(x[1])+whereZero(x[1]-1.)+whereZero(x[2])+whereZero(x[2]-1.)
   print "@ generated %d-dimension mesh with %d elements in each direction"%(d,ne)
   # for ncomp in [1,2]:
   for ncomp in [2]:
        if ncomp==1:
           maskf=1.
           Z=makeRandomFloats(d,-1.,0.)
           K_sup=makeRandomFloatMatrix(d,d,-1.,1.)
           K=numarray.identity(d)*1.
        else:
           maskf=numarray.ones(ncomp)
           Z=makeRandomFloatMatrix(ncomp,d,-1.,0.)
           K_sup=makeRandomFloatTensor(ncomp,d,ncomp,d,-1.,1.)
           K=numarray.zeros([ncomp,d,ncomp,d])*0.
           for i in range(ncomp):
             K[i,:,i,:]=numarray.identity(d)*1.
        K_sup=numarray.array(K_sup)
        Z=numarray.array(Z)
        Z/=length(Z)
        if ncomp==1:
           Zx=Z[0]*x[0]
           for j in range(1,d):
              Zx+=Z[j]*x[j]
        else:
           Zx=x[0]*Z[:,0]
           for j in range(1,d):
              Zx+=x[j]*Z[:,j]
        K+=0.02*K_sup/length(K_sup)
        K/=length(K)
        if ncomp==1:
           U=numarray.matrixmultiply(numarray.transpose(K),Z)
        else:
           U=numarray.zeros([ncomp,d,ncomp])*1.
           for m in range(ncomp):
              for l in range(ncomp):
                 for j in range(d):
                    for i in range(d):
                       U[m,j,l]+=K[m,i,l,j]*Z[l,i]

        # create domain:
        mypde=AdvectivePDE(mydomain)
        # mypde.setSolverMethod(mypde.DIRECT)
        print K
        mypde.setValue(q=msk*maskf)
        mypde.setValue(A=K)
        mypde.setValue(A=K,q=msk*maskf)
        mypde.checkSymmetry()
        # run Peclet
        for Pe in [0.001,1.,1.,10.,100,1000.,10000.,100000.,1000000.,10000000.]:
           peclet=Pe*length(U)/2./length(K)/ne
           print "@@@  Peclet Number :",peclet
           u_ex=exp(Pe*Zx)
           mypde.setValue(r=u_ex)
           # mypde.setValue(B=Data(),C=Pe*U)
           # u=mypde.getSolution()
           # print "@@@@ C=U: Pe = ",peclet,printError(u,u_ex)
           mypde.setValue(C=Data(),B=-Pe*U)
           u=mypde.getSolution()
           print "@@@@ B=-U: Pe = ",peclet,printError(u,u_ex)
