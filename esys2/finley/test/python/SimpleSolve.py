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

print "\nSimpleSolve.py"
print "--------------"

my_options = {
          "verbose" : True,
          "reordering" : NO_REORDERING,
          "tolerance" : 1.E-8,
          "final_residual" : 0.,
          "iterative_method" : PCG,
          "preconditioner" : JACOBI,
          "iter_max" :  5000,
          "iter" : 0,
          "drop_tolerance" : 1.10,
          "drop_storage" : 2.
}

alpha=0.01

print "\nOptions: ", my_options

# generate mesh

print "\nGenerate mesh: finley.Rectangle(9,12,1)=>"
mydomain=finley.Rectangle(9,12,1)

print "\nGenerate mesh: finley.Rectangle(4,4,1)=>"
mydomain=finley.Rectangle(4,4,1)

print "\nGenerate mesh: finley.Rectangle(151,151,1)=>"
mydomain=finley.Rectangle(151,151,1)

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

print "mypde=linearPDE( A=[[1.,0.7],[0.7,1.]], D=alpha, Y=alpha, domain=mydomain)"
mypde=linearPDE(A=[[1.,0.7],[0.7,1.]],D=alpha,Y=alpha,domain=mydomain)
#mypde=linearPDE(A=[[1.,0.],[0.,1.]],D=alpha,Y=alpha,domain=mydomain)

# generate a test solution 1

print "\nGenerate a test solution (1)"
print "----------------------------"

print "\nDirect Solver (1)=>"

u_d=mypde.getSolution(iterative=False,**my_options)

print "\nIterative Solver (1)=>"

u_i=mypde.getSolution(iterative=True,**my_options)

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

print "mypde=linearPDE(A=[[1.,0.],[0.,1.]],q=msk,r=u_ex)"
mypde=linearPDE(A=[[1.,0.],[0.,1.]],q=msk,r=u_ex)

# generate a test solution 2

print "\nGenerate a test solution (2)"
print "----------------------------"

print "\nDirect Solver (2)=>"

u_d=mypde.getSolution(iterative=False,**my_options)

print "\nIterative Solver (2)=>"

u_i=mypde.getSolution(iterative=True,**my_options)

print "\n******************************************************************"
error=u_ex-u_d
print "norm of the error for direct solver is   : ",error.Lsup()/norm_u_ex
error=u_ex-u_i
print "norm of the error for iterative solver is: ",error.Lsup()/norm_u_ex
print "******************************************************************"

print "\n-----"
print "Done."
print "-----"
