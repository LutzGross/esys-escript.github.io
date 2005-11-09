#test_advect.py Gauss


from esys.escript import *
import esys.finley
from esys.linearPDEs import LinearPDE
from esys.linearPDEs import AdvectivePDE
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

if os.path.exists("results/sylvain/dist_error.dat"):
  os.unlink("results/sylvain/dist_error.dat")


xx = mesh.getX()[0]
yy = mesh.getX()[1]

gauss = (160*(xx-0.25)**4 - 320*(xx-0.25)**3 + 160*(xx-0.25)**2 )*( 160*(yy-0.5)**4 - 320*(yy-0.5)**3 + 160*(yy-0.5)**2)
mask_tronc = (xx-0.25).wherePositive()*(1.25-xx).wherePositive()*(yy-0.5).wherePositive()*(1.5-yy).wherePositive()
gauss_tronc_new = gauss*mask_tronc

reference = Lsup(gauss_tronc_new)
gauss_tronc_new.saveDX("results/sylvain/gauss.%2.2i.dx" % 0)

h = Lsup(mesh.getSize())

coeff = 0.3

v_adim = coeff*h/dt
cte = h/(2.0*v_adim)
t_step_end = 2.0/(coeff*h)

velocity = Vector(0.0, Function(mesh))
velocity[0] = v_adim
velocity.saveDX("results/sylvain/velocity_field.dx")


taylor_galerkin_lhs()

q = mesh.getX()[0].whereZero() + (4.0-mesh.getX()[0]).whereZero() + mesh.getX()[1].whereZero() + (2.0-mesh.getX()[1]).whereZero()
surfacePDE.setValue(q=q)


while (t_step <= t_step_end):
   
  gauss_tronc_old = gauss_tronc_new
  
  taylor_galerkin_rhs()
  
  gauss_tronc_new = surfacePDE.getSolution()

  gauss_tronc_old.saveDX("results/sylvain/gauss.%2.2i.dx" % t_step)
  print "integral of f", integrate(gauss_tronc_new*Scalar(1.0, Function(mesh)))
  
  dist = v_adim*dt*t_step
  
  error = 100*abs(Lsup(gauss_tronc_new)-reference)/reference
  File1 = file("results/sylvain/dist_error.dat", "a", 1)
  File1.write("%e %e\n" % (dist, error))
  File1.close()
  t_step = t_step + 1
