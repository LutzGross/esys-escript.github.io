################################################
##                                            ##
## October 2006                               ## 
##                                            ##
##  3D Rayleigh-Taylor instability benchmark  ##
##           by Laurent Bourgouin             ##
##                                            ##
################################################


### IMPORTS ###
from esys.escript import *
import esys.finley
from esys.finley import finley
from esys.escript.linearPDEs import LinearPDE
from esys.escript.pdetools import Projector
import sys
import math

### DEFINITION OF THE DOMAIN ###
l0=1.
l1=1.
n0=20  # IDEALLY 80...
n1=20  # IDEALLY 80...
mesh=esys.finley.Brick(l0=l0, l1=l1, l2=l0, order=2, n0=n0, n1=n1, n2=n0)

### PARAMETERS OF THE SIMULATION ###
rho1 = 1.0e3         # DENSITY OF THE FLUID AT THE BOTTOM
rho2 = 1.01e3        # DENSITY OF THE FLUID ON TOP
eta1 = 1.0e2         # VISCOSITY OF THE FLUID AT THE BOTTOM
eta2 = 1.0e2         # VISCOSITY OF THE FLUID ON TOP
penalty = 1.0e3      # PENALTY FACTOR FOT THE PENALTY METHOD
g=10.                # GRAVITY
t_step = 0
t_step_end = 2000
reinit_max = 30      # NUMBER OF ITERATIONS DURING THE REINITIALISATION PROCEDURE
reinit_each = 3      # NUMBER OF TIME STEPS BETWEEN TWO REINITIALISATIONS
h = Lsup(mesh.getSize())
numDim = mesh.getDim()
smooth = h*2.0       # SMOOTHING PARAMETER FOR THE TRANSITION ACROSS THE INTERFACE

### DEFINITION OF THE PDE ###
velocityPDE = LinearPDE(mesh, numEquations=3)
velocityPDE.setSolverMethod(solver=LinearPDE.DIRECT)
advectPDE = LinearPDE(mesh)
advectPDE.setReducedOrderOn()
advectPDE.setValue(D=1.0)
advectPDE.setSolverMethod(solver=LinearPDE.DIRECT)
reinitPDE = LinearPDE(mesh, numEquations=1)
reinitPDE.setReducedOrderOn()
reinitPDE.setSolverMethod(solver=LinearPDE.LUMPING)
my_proj=Projector(mesh)

### BOUNDARY CONDITIONS ###
xx = mesh.getX()[0]
yy = mesh.getX()[1]
zz = mesh.getX()[2]
top = whereZero(zz-l1)
bottom = whereZero(zz)
left = whereZero(xx)
right = whereZero(xx-l0)
front = whereZero(yy)
back = whereZero(yy-l0)
b_c = (bottom+top)*[1.0, 1.0, 1.0] + (left+right)*[1.0,0.0, 0.0] + (front+back)*[0.0, 1.0, 0.0]
velocityPDE.setValue(q = b_c)

pressure = Scalar(0.0, ContinuousFunction(mesh))

### INITIALISATION OF THE INTERFACE ###
func = -(-0.1*cos(math.pi*xx/l0)*cos(math.pi*yy/l0)-zz+0.4)
phi = func.interpolate(ReducedSolution(mesh))


def advect(phi, velocity, dt):
### SOLVES THE ADVECTION EQUATION ###
 
  Y = phi.interpolate(Function(mesh))
  for i in range(numDim):
    Y -= (dt/2.0)*velocity[i]*grad(phi)[i]
  advectPDE.setValue(Y=Y)    
  phi_half = advectPDE.getSolution()

  Y = phi
  for i in range(numDim):
    Y -= dt*velocity[i]*grad(phi_half)[i]
  advectPDE.setValue(Y=Y)    
  phi = advectPDE.getSolution()

  print "Advection step done"
  return phi

def reinitialise(phi):
### SOLVES THE REINITIALISATION EQUATION ###
  s = sign(phi.interpolate(Function(mesh)))
  w = s*grad(phi)/length(grad(phi))
  dtau = 0.3*h
  iter =0
  previous = 100.0
  mask = whereNegative(abs(phi)-1.2*h)
  reinitPDE.setValue(q=mask, r=phi)
  while (iter<=reinit_max):
    prod_scal =0.0
    for i in range(numDim):
      prod_scal += w[i]*grad(phi)[i]
    coeff = s - prod_scal
    ps2=0
    for i in range(numDim):
      ps2 += w[i]*grad(my_proj(coeff))[i]
    reinitPDE.setValue(D=1.0, Y=phi+dtau*coeff-0.5*dtau**2*ps2)
    phi = reinitPDE.getSolution()
    error = Lsup((previous-phi)*whereNegative(abs(phi)-3.0*h))/h
    print "Reinitialisation iteration :", iter, " error:", error
    previous = phi
    iter +=1
  return phi

def update_phi(phi, velocity, dt, t_step):
### CALLS THE ADVECTION PROCEDURE AND THE REINITIALISATION IF NECESSARY ###  
  phi=advect(phi, velocity, dt)
  if t_step%reinit_each ==0:
    phi = reinitialise(phi)
  return phi
  
def update_parameter(phi, param_neg, param_pos):
### UPDATES THE PARAMETERS TABLE USING THE SIGN OF PHI, A SMOOTH TRANSITION IS DONE ACROSS THE INTERFACE ###
  mask_neg = whereNonNegative(-phi-smooth)
  mask_pos = whereNonNegative(phi-smooth)
  mask_interface = whereNegative(abs(phi)-smooth)
  param = param_pos*mask_pos + param_neg*mask_neg + ((param_pos+param_neg)/2 +(param_pos-param_neg)*phi/(2.*smooth))*mask_interface
  return param

def solve_vel(rho, eta, pressure):
### SOLVES THE VELOCITY PROBLEM USING A PENALTY METHOD FOR THE INCOMPRESSIBILITY ###
  error = 1.0
  ref = pressure*1.0
  p_iter=0
  while (error >= 1.0e-2):
  
    A=Tensor4(0.0, Function(mesh))
    for i in  range(numDim):
      for j in range(numDim):
        A[i,j,i,j] += eta
        A[i,j,j,i] += eta
	A[i,i,j,j] += penalty*eta

    Y = Vector(0.0,Function(mesh))
    Y[1] -= rho*g
    
    X = Tensor(0.0, Function(mesh))
    for i in range(numDim):
      X[i,i] += pressure
    
    velocityPDE.setValue(A=A, X=X, Y=Y)
    velocity_new = velocityPDE.getSolution()
    p_iter +=1
    if p_iter >=500:
      print "You're screwed..."
      sys.exit(1)    
    
    pressure -= penalty*eta*(trace(grad(velocity_new)))
    error = penalty*Lsup(trace(grad(velocity_new)))/Lsup(grad(velocity_new))
    print "\nPressure iteration number:", p_iter
    print "error", error
    ref = pressure*1.0
    
  return velocity_new, pressure
  
### MAIN LOOP, OVER TIME ###
while t_step <= t_step_end:

  rho = update_parameter(phi, rho1, rho2)
  eta = update_parameter(phi, eta1, eta2)

  velocity_new, pressure = solve_vel(rho, eta, pressure)
  dt = 0.3*Lsup(mesh.getSize())/Lsup(velocity_new)
  phi = update_phi(phi, velocity_new, dt, t_step)

### PSEUDO POST-PROCESSING ###
  print "##########  Saving image", t_step, " ###########" 
  phi.saveVTK("/home/laurent/results2006/instability/phi3D.%2.2i.vtk" % t_step)  

  print "######################"
  print "Time step:", t_step
  print "######################"
  t_step += 1
