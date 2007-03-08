# $Id$

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

print "\nSimpleSolve.py"
print "--------------"

alpha=0.7
error_tol=1.e-5

# generate mesh

# print "\nGenerate mesh: finley.Rectangle(9,12,1)=>"
# mydomain=finley.Rectangle(140,140)

# print "\nGenerate mesh: finley.Rectangle(4,4,1)=>"
# mydomain=finley.Rectangle(5,5,1)
print "A"
mydomain=finley.Rectangle(500,500,1)
print "B"
mydomain=finley.Rectangle(150,150,1)
print "C"

print "\nGenerate mesh: finley.Rectangle(151,151,1)=>"
# mydomain=finley.Rectangle(151,151,1)
# mydomain=finley.Rectangle(128,128,1)

# set the direct solver switch
DIRECT=LinearPDE.DIRECT
# DIRECT=LinearPDE.ITERATIVE # this will switch of the DIRECT solver call to avoid external library calls which may not be available everywhere

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
print "is x protected: ",x.isProtected()


print "norm_u_ex=Lsup(u_ex):"
norm_u_ex=Lsup(u_ex)
print "is u_ex protected: ",u_ex.isProtected()

print "\nGenerate a test solution (1)"
print "----------------------------"

print "mypde=LinearPDE( A=[[1.,0.8],[0.4,1.]], D=alpha, Y=alpha, domain=mydomain)"
mypde=LinearPDE(mydomain)
mypde.setDebugOn()
mypde.setValue(A=[[1.,-0.001],[-0.001,1.]],D=alpha,Y=alpha)

print "mypde.checkSymmetry()"
print mypde.checkSymmetry()

print "\nIterative Solver (1)=>"
# mypde.setSolverMethod(mypde.PRES20,preconditioner=mypde.ILU0)
mypde.setSolverMethod(mypde.BICGSTAB,preconditioner=mypde.JACOBI)
u_i=mypde.getSolution(verbose=True,iter_max=3000)

print "\nDirect Solver (1)=>"
mypde.setSolverMethod(DIRECT)
u_d=mypde.getSolution(verbose=True)

print "\n***************************************************************"
error=u_ex-u_d
error_norm=Lsup(error)/norm_u_ex
print "norm of the error for direct solver is   : ",error_norm
if error_norm > error_tol:
  print "### error norm exceeded maximum tolerance ###"
  sys.exit(1)
error=u_ex-u_i
error_norm=Lsup(error)/norm_u_ex
print "norm of the error for iterative solver is: ",error_norm
if error_norm > error_tol:
  print "### error norm exceeded maximum tolerance ###"
  sys.exit(1)
print "***************************************************************"
# del mypde
del u_i
# releaseUnusedMemory()
print "***************************************************************"



# get handles to nodes and elements 2

print "\nGet handles to nodes and elements(2)=>"
print "--------------------------------------"

print "x=n.getX():"
x=n.getX()

print " msk=whereZero(x[0])+whereZero(x[0]-1.)"
msk=whereZero(x[0])+whereZero(x[0]-1.)

print "mypde=LinearPDE(A=[[1.,0.],[0.,1.]],q=msk,r=u_ex)"
mypde=LinearPDE(mydomain)
mypde.setDebugOn()
mypde.setValue(A=[[1.,0.0],[0.0,1.]],q=msk,r=u_ex)

print "mypde.checkSymmetry()"
print mypde.checkSymmetry()

# generate a test solution 2

print "\nGenerate a test solution (2)"
print "----------------------------"

print "\nDirect Solver (2)=>"

mypde.setSymmetryOn() 
mypde.setTolerance(1.e-13)

# mypde.setSymmetryOn() : is not woking yet!
mypde.setSolverMethod(DIRECT)
u_d=mypde.getSolution(verbose=True)

print "\nIterative Solver (2)=>"

mypde.setSolverMethod(mypde.PCG)
u_i=mypde.getSolution(verbose=True,iter_max=3000)

print "\n******************************************************************"
error=u_ex-u_d
error_norm=Lsup(error)/norm_u_ex
print "norm of the error for direct solver is   : ",error_norm
if error_norm > error_tol:
  print "### error norm exceeded maximum tolerance ###"
  sys.exit(1)
error=u_ex-u_i
error_norm=Lsup(error)/norm_u_ex
print "norm of the error for iterative solver is: ",error_norm
if error_norm >  error_tol:
  print "### error norm exceeded maximum tolerance ###"
  sys.exit(1)
print "******************************************************************"

print "\n-----"
print "Done."
print "-----"

stoptime = time.time()
elapsed = stoptime - starttime
print "\nElapsed time: ", elapsed, "\n"

sys.exit(0)
