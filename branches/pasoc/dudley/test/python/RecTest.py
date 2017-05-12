
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

import sys
import os
import esys.escriptcore.utestselect as unittest

from esys.escript import *
from esys.escript.linearPDEs import *
from esys import dudley

import numpy

Pi=numpy.pi
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
    mypde.setSymmetryOn()
    mypde.setDebugOn()
    mypde.setReducedOrderTo(reduce)
    mypde.setValue(A=numpy.identity(msh.getDim()),D=sml,Y=(sml+4*Pi**2*msh.getDim())*u_ex,q=constraints,r=u_ex)
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
              msh=dudley.Rectangle(numElements,numElements,order,periodic0=i0,periodic1=i1,useElementsOnFace=onElements)
              n=ContinuousFunction(msh)
              x=n.getX()
              c=Scalar(0,what=n) 
              if i0==False:
                  c+=whereZero(x[0])+whereZero(x[0]-1.)
              if i1==False:
                  c+=whereZero(x[1])+whereZero(x[1]-1.)
              error=TheTest(msh,c,reduce)
              text="Rectangle order = %d%s%s, periodic0= %d, periodic1= %d: error= %f"%(order,redtext,onElmtext,i0,i1,error)
              print("@@ ",text)
              if error>max_error: 
                  max_error=error
                  max_text=text
        elif dim==3:
          for i0 in [True,False]:
            for i1 in [True,False]:
              for i2 in [True,False]:
                msh=dudley.Brick(numElements,numElements,numElements,order,periodic0=i0,periodic1=i1,periodic2=i2,useElementsOnFace=onElements)
                n=ContinuousFunction(msh)
                x=n.getX()
                c=Scalar(0,what=n)
                if i0==False:
                  c+=whereZero(x[0])+whereZero(x[0]-1.)
                if i1==False:
                  c+=whereZero(x[1])+whereZero(x[1]-1.)
                if i2==False:
                  c+=whereZero(x[2])+whereZero(x[2]-1.)
                error=TheTest(msh,c,reduce)
                text="Brick order = %d%s%s, periodic0= %d, periodic1= %d, periodic2= %d: error= %f"%(order,redtext,onElmtext,i0,i1,i2,error)
                print("@@ ",text)
                if error>max_error: 
                    max_error=error
                    max_text=text

print("@@@@ maximum error for :",max_text)
