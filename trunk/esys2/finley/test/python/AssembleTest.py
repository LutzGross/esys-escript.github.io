# assemble test: 
#
# $Id$

import sys
import os
import unittest
                                                                                                                                                          
esys_root=os.getenv('ESYS_ROOT')
sys.path.append(esys_root+'/finley/lib')
sys.path.append(esys_root+'/escript/lib')
sys.path.append(esys_root+'/escript/py_src')
                                                                                                                                                          
from escript import *
from util import *
from linearPDE import *

import finley
from math import *
global seed
num_elem=2 # number of elements in each spatial direction
num_equations=3 # number of equations
integration_order=-1 # order of the integration scheme
order_u=1 # solution order


seed=1
rand_factor=sqrt(2.)
#
#   The test solution is represented in the basis 
#
#   B=[1,x_1,x_2,x_3,x_1**2,x_2**2,x_3**2,....,x_1**order_u,x_2**order_u,x_3**order_u] for dim=3
#                 or 
#   B=[1,x_1,x_2,x_1**2,x_2**2,....,x_1**order_u,x_2**order_u]   for dim=2
#
#   length of B is len_Basis.
#
#   any test solution u with numComp components is represented by the numComp x len_Basis matrix sol.
#   the actual solution is then given by u=matmul(sol,B). 
#

global maxError,coeff,total_maxError,total_coeff
total_maxError=0
total_coeff=""
maxError=0
coeff=""
D=()


def randomNum():
     global seed
     seed+=1
     s=rand_factor*seed
     return s-int(s)
    
def generateRandom(arg):
    if len(arg)==0:
       return randomNum()
    out=numarray.ones(arg,numarray.Float64) # *randomNum()
    return out

def algebraicGrad(u):
   out=numarray.zeros([u.shape[0],dim,len_Basis],numarray.Float64)
   for s in range(u.shape[0]):
     for i in range(dim):
        for k in range(len_Basis):
          h=0
          for j in range(len_Basis):
             h+=u[s,j]*D[i][k,j]
          out[s,i,k]=h
   return out

def algebraicDiv(u):
   out=numarray.zeros([u.shape[0],len_Basis],numarray.Float64)
   for s in range(u.shape[0]):
     for k in range(len_Basis):
        h=0
        for i in range(dim):
          for j in range(len_Basis):
             h+=u[s,i,j]*D[i][k,j]
        out[s,k]=h
   return out

def mult4(A,u):
      out=numarray.zeros([A.shape[0],A.shape[1],len_Basis],numarray.Float64)
      for s in range(A.shape[0]):
        for i in range(A.shape[1]):
           for k in range(len_Basis):
             h=0
             for t in range(A.shape[2]):
               for j in range(A.shape[3]):
                 h+=A[s,i,t,j]*u[t,j,k]
             out[s,i,k]=h
      return out

def mult3_2(A,u):
      out=numarray.zeros([A.shape[0],A.shape[1],len_Basis],numarray.Float64)
      for s in range(A.shape[0]):
        for i in range(A.shape[1]):
           for k in range(len_Basis):
             h=0
             for t in range(A.shape[2]):
                 h+=A[s,i,t]*u[t,k]
             out[s,i,k]=h
      return out

def mult3_1(A,u):
      out=numarray.zeros([A.shape[0],len_Basis],numarray.Float64)
      for s in range(A.shape[0]):
           for k in range(len_Basis):
             h=0
             for t in range(A.shape[1]):
               for j in range(A.shape[2]):
                 h+=A[s,t,j]*u[t,j,k]
             out[s,k]=h
      return out

def mult2(A,u):
      out=numarray.zeros([A.shape[0],len_Basis],numarray.Float64)
      for s in range(A.shape[0]):
         for k in range(len_Basis):
           h=0
           for t in range(A.shape[1]):
               h+=A[s,t]*u[t,k]
           out[s,k]=h
      return out


def eval(u,this):
    x=this.getX()
    if u.rank==2:
       out=Data(value=0,shape=(u.shape[0],),what=this,expand=True)
       for i0 in range(u.shape[0]):
          out[i0]=u[i0,0]
          for p in range(order_u):
             for j in range(dim):
               out[i0]+=u[i0,p*dim+j+1]*x[j]**(p+1)
    else:
       if u.shape[0]==1:
          out=Data(value=0,shape=(u.shape[1],),what=this,expand=True)
          for i1 in range(u.shape[1]):
             out[i1]=u[0,i1,0]
             for p in range(order_u):
                for j in range(dim):
                  out[i1]+=u[0,i1,p*dim+j+1]*x[j]**(p+1)

       elif u.shape[1]==1:
          out=Data(value=0,shape=(u.shape[0],),what=this,expand=True)
          for i0 in range(u.shape[0]):
             out[i0]=u[i0,0,0]
             for p in range(order_u):
                for j in range(dim):
                  out[i0]+=u[i0,0,p*dim+j+1]*x[j]**(p+1)
       else:
          out=Data(value=0,shape=(u.shape[0],u.shape[1]),what=this,expand=True)
          for i0 in range(u.shape[0]):
           for i1 in range(u.shape[1]):
             out[i0,i1]=u[i0,i1,0]
             for p in range(order_u):
                for j in range(dim):
                  out[i0,i1]+=u[i0,i1,p*dim+j+1]*x[j]**(p+1)
    return out

def checkSystem(text,operator,u,rhs):
    global maxError,coeff
    error=Lsup(operator*u-rhs)
    print "@@ "+text+" error: ",error
    if error>=maxError:
       maxError=error
       coeff=text
#
#
def TestSystem(numEqu,numComp,mydomain,reduce):
   elem=Function(mydomain)
   face_elem=FunctionOnBoundary(mydomain)
   nodes=ContinuousFunction(mydomain)
   nrml=face_elem.getNormal()
   #
   #   test solution:
   #
   u=generateRandom([numComp,len_Basis])
   # u=numarray.zeros([numComp,len_Basis])
   # u[0,0]=1
   # u[0,1]=1
   U=eval(u,nodes)
   gradu=algebraicGrad(u)
   #
   #   test A:
   #
   for p in range(numEqu):
      for q in range(numComp):
         for i in range(dim):
            for j in range(dim):

               # check div( A grad(u) ) = div( X )
               c_A=numarray.zeros([numEqu,dim,numComp,dim])
               c_A[p,i,q,j]=1
               if numEqu==1 and numComp==1:
                  c_A2=numarray.reshape(c_A,[dim,dim])
                  text="A[%d,%d]"%(i,j)
               else:
                  c_A2=c_A
                  text="A[%d,%d,%d,%d]"%(p,i,q,j)
               x=mult4(c_A,gradu)
               mypde1=linearPDE(domain=mydomain,A=c_A2,X=eval(x,elem))
               mypde1.setReducedOrderForSolutionsTo(reduce)
               checkSystem(text+" const with X",mypde1.getOperator(),U,mypde1.getRightHandSide())
               # check div( A grad(u) ) = Y
               y=-algebraicDiv(x)
               mypde2=linearPDE(domain=mydomain,Y=eval(y,elem),y=matmult(eval(x,face_elem),nrml))
               mypde2.setReducedOrderForSolutionsTo(reduce)
               checkSystem(text+" const with Y",mypde1.getOperator(),U,mypde2.getRightHandSide())

            # check div( B u ) = div( X )
            c_B=numarray.zeros([numEqu,dim,numComp])
            c_B[p,i,q]=1
            if numEqu==1 and numComp==1:
               c_B2=numarray.reshape(c_B,[dim])
               text="B[%d]"%(i)
            else:
               c_B2=c_B
               text="B[%d,%d,%d]"%(p,i,q)
            x=mult3_2(c_B,u)
            mypde1=linearPDE(domain=mydomain,B=c_B2,X=eval(x,elem))
            mypde1.setReducedOrderForSolutionsTo(reduce)
            checkSystem(text+" const with X",mypde1.getOperator(),U,mypde1.getRightHandSide())
            # check div( B u ) = Y
            y=-algebraicDiv(x)
            mypde2=linearPDE(domain=mydomain,Y=eval(y,elem),y=matmult(eval(x,face_elem),nrml))
            mypde2.setReducedOrderForSolutionsTo(reduce)
            checkSystem(text+" const with Y",mypde1.getOperator(),U,mypde2.getRightHandSide())

            # check C grad(u) = Y
            c_C=numarray.zeros([numEqu,numComp,dim])
            c_C[p,q,i]=1
            if numEqu==1 and numComp==1:
               c_C2=numarray.reshape(c_C,[dim])
               text="C[%d]"%(i)
            else:
               c_C2=c_C
               text="C[%d,%d,%d]"%(p,q,i)
            y=mult3_1(c_C,gradu)
            mypde1=linearPDE(domain=mydomain,C=c_C2,Y=eval(y,elem))
            mypde1.setReducedOrderForSolutionsTo(reduce)
            checkSystem(text+" const with Y",mypde1.getOperator(),U,mypde1.getRightHandSide())


         # check D u= Y
         c_D=numarray.zeros([numEqu,numComp])
         c_D[p,q]=1
         if numEqu==1 and numComp==1:
            c_D2=numarray.reshape(c_D,[1])
            text="D"
         else:
            c_D2=c_D
            text="D[%d,%d]"%(p,q)
         y=mult2(c_D,u)
         mypde1=linearPDE(domain=mydomain,D=c_D2,Y=eval(y,elem))
         mypde1.setReducedOrderForSolutionsTo(reduce)
         checkSystem(text+" const with Y",mypde1.getOperator(),U,mypde1.getRightHandSide())

#  we start:
for dim in [2,3]:
    len_Basis=dim*order_u+1
    #
    #  the differential operators:
    #   
    #  D_j(matmul(sol,B))=matmul(sol,D_jB)=matmul(matmul(sol,D[j]),B)
    #
    D=()
    for i in range(dim):
      D=D+(numarray.zeros([len_Basis,len_Basis],numarray.Float64),)
      for j in range(order_u):
         if j==0:
           D[i][0,i+j+1]=1
         else:
           D[i][(j-1)*dim+i+1,j*dim+i+1]=j+1
    #
    #   generate mydomain:
    #
    for order in [1,2]:
     for onElements in [False,True]:
        if onElements==True:
            onElmtext=", with elements on faces"
        else:
            onElmtext=""
        if dim==3:
           mydomain=finley.Brick(num_elem,num_elem,num_elem,order=order,integrationOrder=integration_order,useElementsOnFace=onElements)
        else:
           mydomain=finley.Rectangle(num_elem,num_elem,order=order,integrationOrder=integration_order,useElementsOnFace=onElements)
        for reduce in [False,True]:
           if reduce==True:
              redtext=",reduced"
           else:
              redtext=""
           #
           #  and start the test process:
           #
           for numEqu in range(1,num_equations+1):
               print "@@@ Start testing assembling with dim=%d and %d equations (order=%d%s%s)"%(dim,numEqu,order,redtext,onElmtext)
               TestSystem(numEqu,numEqu,mydomain,reduce)
               print "@@@ end testing assembling (order=%d%s%s) with %d equations in %dD with maximum error %e for %s"%(order,redtext,onElmtext,numEqu,dim,maxError,coeff)
               if maxError>=total_maxError:
                  total_maxError=maxError
                  total_coeff=coeff+" with %d equations in %dD (order=%d%s%s)"%(numEqu,dim,order,redtext,onElmtext)
    
print "@@@@ end testing assemblage with maximal error %e for %s"%(total_maxError,total_coeff)


