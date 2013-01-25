
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
__url__="https://launchpad.net/escript-finley"

#test_advect.py Gauss

from esys.escript import *
import esys.finley
from esys.escript.linearPDEs import LinearPDE
from esys.escript.linearPDEs import AdvectivePDE
import os

def classical_lhs():
#######################   FINLEY UPWINDING HARDCODED    #################
  global surfacePDE
  surfacePDE = LinearPDE(mesh)
  
  A = Tensor(0.0, Function(mesh))
  for j in range(2):
    for l in range(2):
      A[j,l] = velocity[j] * velocity[l] * cte * dt
     
  B = velocity * cte
  C = dt * velocity
  D=1.0
  
  surfacePDE.setValue(A=A, B=B, C=C,  D=D)

def classical_rhs():
#######################   FINLEY UPWINDING HARDCODED    #################
  X = velocity * cte * gauss_tronc_old
  Y = gauss_tronc_old
  surfacePDE.setValue(X=X, Y=Y)
  
def finley_upwd_lhs():
#######################   FINLEY UPWINDING   #################
  global surfacePDE
  surfacePDE = AdvectivePDE(mesh)
  
  C = dt * velocity
  D=1.0
  surfacePDE.setValue(C=C,  D=D)
  
def finley_upwd_rhs():
#######################   FINLEY UPWINDING   #################
  Y = gauss_tronc_old
  surfacePDE.setValue(Y=Y)

def taylor_galerkin_lhs():
#######################   TAYLOR - GALERKIN   #################
  global surfacePDE
  surfacePDE = LinearPDE(mesh)
  
  surfacePDE.setValue(D=1.0)
  
def taylor_galerkin_rhs():
#######################   TAYLOR - GALERKIN   #################
  X = Vector(0.0, Function(mesh))
  for i in range(2):
    for j in range(2):
      X[i] -= (dt**2)/2.0*velocity[i]*velocity[j]*grad(gauss_tronc_old)[j]
  Y = gauss_tronc_old*Scalar(1.0, Function(mesh))
  for j in range(2):
    Y -= dt*velocity[j]*grad(gauss_tronc_old)[j]
    
  surfacePDE.setValue(X=X, Y=Y)


############ start of main code ########################


dt = 0.5
t_step = 1


mesh = esys.finley.Rectangle(l0=4.0, l1=2.0, order=1, n0=60, n1=30)


xx = mesh.getX()[0]
yy = mesh.getX()[1]

gauss = (160*(xx-0.25)**4 - 320*(xx-0.25)**3 + 160*(xx-0.25)**2 )*( 160*(yy-0.5)**4 - 320*(yy-0.5)**3 + 160*(yy-0.5)**2)
mask_tronc = wherePositive(xx-0.25)*wherePositive(1.25-xx)*wherePositive(yy-0.5)*wherePositive(1.5-yy)
gauss_tronc_new = gauss*mask_tronc

reference = Lsup(gauss_tronc_new)

h = Lsup(mesh.getSize())

coeff = 0.3

v_adim = coeff*h/dt
cte = h/(2.0*v_adim)
t_step_end = 2.0/(coeff*h)

velocity = Vector(0.0, Function(mesh))
velocity[0] = v_adim


taylor_galerkin_lhs()

q = whereZero(mesh.getX()[0]) + whereZero(4.0-mesh.getX()[0]) + whereZero(mesh.getX()[1]) + whereZero(2.0-mesh.getX()[1])
surfacePDE.setValue(q=q)


while (t_step <= t_step_end):
   
  gauss_tronc_old = gauss_tronc_new
  
  taylor_galerkin_rhs()
  
  gauss_tronc_new = surfacePDE.getSolution()

  print "integral of f", integrate(gauss_tronc_new*Scalar(1.0, Function(mesh)))
  
  dist = v_adim*dt*t_step
  
  error = 100*abs(Lsup(gauss_tronc_new)-reference)/reference
  t_step = t_step + 1
