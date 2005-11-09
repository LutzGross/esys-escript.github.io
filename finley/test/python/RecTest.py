# $Id$

import sys
import os
import unittest

from esys.escript import *
from esys.escript.linearPDEs import *
from esys import finley

import numarray

Pi=numarray.pi
numElements=15
sml=0.1

#
# this test the assemblage of problems with periodic boundary conditions:
# 

#
#    test solution is u_ex=sin(2*Pi*n*x0)*...*sin(2*Pi*n*x_dim)
#
def TheTest(msh,constraints,reduce):
    n=ContinuousFunction(msh)
    x=n.getX()
    # set a test problem:
    u_ex=Scalar(1,what=n)
    for i in range(msh.getDim()):
      u_ex*=sin(2*Pi*x[i])
    mypde=LinearPDE(msh)
    mypde.setValue(A=numarray.identity(msh.getDim()),D=sml,Y=(sml+4*Pi**2*msh.getDim())*u_ex,q=constraints,r=u_ex)
    mypde.setSymmetryOn()
    mypde.setDebugOn()
    mypde.setReducedOrderTo(reduce)
    return Lsup(mypde.getSolution()-u_ex)/Lsup(u_ex)

max_error=0
max_text="none"
for onElements in [False,True]:
 if onElements==True:
    onElmtext=", with elements on faces"
 else:
    onElmtext=""
 for order in [1,2]:
  for reduce in [False,True]:
    if reduce==True:
       redtext=",reduced"
    else:
       redtext=""
    for dim in [2,3]:
        if dim==2:
          for i0 in [True,True]:
            for i1 in [True,True]:
              msh=finley.Rectangle(numElements,numElements,order,periodic0=i0,periodic1=i1,useElementsOnFace=onElements)
              n=ContinuousFunction(msh)
              x=n.getX()
              c=Scalar(0,what=n) 
              if i0==False:
                  c+=x[0].whereZero()+(x[0]-1.).whereZero()
              if i1==False:
                  c+=x[1].whereZero()+(x[1]-1.).whereZero()
              error=TheTest(msh,c,reduce)
              text="Rectangle order = %d%s%s, periodic0= %d, periodic1= %d: error= %f"%(order,redtext,onElmtext,i0,i1,error)
              print "@@ ",text
              if error>max_error: 
                  max_error=error
                  max_text=text
        elif dim==3:
          for i0 in [True,False]:
            for i1 in [True,False]:
              for i2 in [True,False]:
                msh=finley.Brick(numElements,numElements,numElements,order,periodic0=i0,periodic1=i1,periodic2=i2,useElementsOnFace=onElements)
                n=ContinuousFunction(msh)
                x=n.getX()
                c=Scalar(0,what=n)
                if i0==False:
                  c+=x[0].whereZero()+(x[0]-1.).whereZero()
                if i1==False:
                  c+=x[1].whereZero()+(x[1]-1.).whereZero()
                if i2==False:
                  c+=x[2].whereZero()+(x[2]-1.).whereZero()
                error=TheTest(msh,c,reduce)
                text="Brick order = %d%s%s, periodic0= %d, periodic1= %d, periodic2= %d: error= %f"%(order,redtext,onElmtext,i0,i1,i2,error)
                print "@@ ",text
                if error>max_error: 
                    max_error=error
                    max_text=text

print "@@@@ maximum error for :",max_text
