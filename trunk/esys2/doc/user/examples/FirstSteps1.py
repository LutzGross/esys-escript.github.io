# import ESyS and finley
from ESyS import *
import finley
# set a value for alpha:
alpha=10
# generate mesh:
mydomain=finley.Rectangle(n0=40,n1=20,l0=2.,l1=1.)
# generate a system:
mypde=linearPDE(A=[[1,0],[0,1]],D=alpha,Y=10,domain=mydomain)
# generate a test solution:
u=mypde.getSolution()
# calculate the error of the solution
error=u-1.
print "norm of the approximation error is ",Lsup(error)
