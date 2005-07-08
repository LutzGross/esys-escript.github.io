# $Id$

import sys
import os
import unittest

from esys.escript import *
from esys.linearPDEs import *
from esys import finley

print "\nSimpleSolve.py"
print "--------------"

alpha=0.025

# generate mesh

# print "\nGenerate mesh: finley.Rectangle(9,12,1)=>"
# mydomain=finley.Rectangle(140,140)

# print "\nGenerate mesh: finley.Rectangle(4,4,1)=>"
# mydomain=finley.Rectangle(5,5,1)

print "\nGenerate mesh: finley.Rectangle(151,151,1)=>"
mydomain=finley.Rectangle(151,151,1)
# mydomain=finley.Rectangle(128,128,1)

print "\nSetup domain and functions"
print "--------------------------"

print "e=Function(mydomain):"
e=Function(mydomain)

print "n=ContinuousFunction(mydomain):"
n=ContinuousFunction(mydomain)

# get handles to nodes and elements 1

print "\nGet handles to nodes and elements(1)=>"
print "--------------------------------------"

print "u_ex=Scalar(1,n,True):"
u_ex=Scalar(1,n,True)

print "x=e.getX():"
x=e.getX()

print "norm_u_ex=u_ex.Lsup():"
norm_u_ex=u_ex.Lsup()

print "\nGenerate a test solution (1)"
print "----------------------------"

print "mypde=LinearPDE( A=[[1.,0.8],[0.4,1.]], D=alpha, Y=alpha, domain=mydomain)"
mypde=LinearPDE(mydomain)
mypde.setDebugOn()
mypde.setValue(A=[[1.,0.8],[0.4,1.]],D=alpha,Y=alpha)

print "mypde.checkSymmetry()"
print mypde.checkSymmetry()

print "\nIterative Solver (1)=>"
# u_i=mypde.getSolution(preconditioner=ILU0,iter_max=3000)
u_i=mypde.getSolution(iter_max=3000)


print "\nDirect Solver (1)=>"
mypde.setSolverMethod(DIRECT)
u_d=mypde.getSolution()

print "\n***************************************************************"
error=u_ex-u_d
print "norm of the error for direct solver is   : ",error.Lsup()/norm_u_ex
error=u_ex-u_i
print "norm of the error for iterative solver is: ",error.Lsup()/norm_u_ex
print "***************************************************************"

# get handles to nodes and elements 2

print "\nGet handles to nodes and elements(2)=>"
print "--------------------------------------"

print "x=n.getX():"
x=n.getX()

print "msk=x[0].whereZero()+(x[0]-1.).whereZero()"
msk=x[0].whereZero()+(x[0]-1.).whereZero()

print "mypde=LinearPDE(A=[[1.,0.],[0.,1.]],q=msk,r=u_ex)"
mypde=LinearPDE(mydomain)
mypde.setDebugOn()
mypde.setValue(A=[[1.,0.],[0.,1.]],q=msk,r=u_ex)

print "mypde.checkSymmetry()"
print mypde.checkSymmetry()

# generate a test solution 2

print "\nGenerate a test solution (2)"
print "----------------------------"

print "\nDirect Solver (2)=>"

#mypde.setSymmetryOn() 
mypde.setTolerance(1.e-13)

# mypde.setSymmetryOn() : is not woking yet!
mypde.setSolverMethod(mypde.DIRECT)
u_d=mypde.getSolution()

print "\nIterative Solver (2)=>"

mypde.setSymmetryOn() 
mypde.setSolverMethod(mypde.DEFAULT_METHOD)
u_i=mypde.getSolution(iter_max=3000)

print "\n******************************************************************"
error=u_ex-u_d
print "norm of the error for direct solver is   : ",error.Lsup()/norm_u_ex
error=u_ex-u_i
print "norm of the error for iterative solver is: ",error.Lsup()/norm_u_ex
print "******************************************************************"

print "\n-----"
print "Done."
print "-----"
