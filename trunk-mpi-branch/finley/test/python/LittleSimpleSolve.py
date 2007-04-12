# $Id: SimpleSolve.py 1016 2007-03-08 06:31:28Z gross $

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import sys
import os
import unittest
import time

from esys.escript import *
from esys.escript.linearPDEs import *
from esys import finley

starttime = time.time()

error_tol=1.e-5
# generate mesh
print "\nGenerate mesh: finley.Rectangle()=>"
# mydomain=finley.Rectangle(150,150,1)
mydomain=finley.Rectangle(80,10,1)

print "e=Function(mydomain):"
e=Function(mydomain)

print "n=ContinuousFunction(mydomain):"
n=ContinuousFunction(mydomain)

print "u_ex=Scalar(1,n,True):"
u_ex=Scalar(1,n,True)

print "x=e.getX():"
x=e.getX()
print "is x protected ? ",x.isProtected()


norm_u_ex=Lsup(u_ex)
print "norm_u_ex=Lsup(u_ex):", norm_u_ex
print "is u_ex protected: ",u_ex.isProtected()

print "x=n.getX():"
x=n.getX()

print " msk=whereZero(x[0])+whereZero(x[0]-1.)"
msk=whereZero(x[0])+whereZero(x[0]-1.)

print "mypde=LinearPDE(A=[[1.,0.],[0.,1.]],q=msk,r=u_ex)"
mypde=LinearPDE(mydomain)
mypde.setDebugOn()
mypde.setValue(A=[[1.,0.0],[0.0,1.]],q=msk,r=u_ex)
mypde.setSymmetryOn()
mypde.setTolerance(1.e-13)
print "enter iterative solver"
mypde.setSolverMethod(mypde.PCG)
u_i=mypde.getSolution(verbose=True,iter_max=3000)
print "leave iterative solver"

error=u_ex-u_i
error_norm=Lsup(error)/norm_u_ex
print "norm of the error for iterative solver is: ",error_norm
if error_norm >  error_tol:
  print "### error norm exceeded maximum tolerance ###"
  sys.exit(1)

stoptime = time.time()
elapsed = stoptime - starttime
print "\nElapsed time: ", elapsed, "\n"

# sys.exit(0)
