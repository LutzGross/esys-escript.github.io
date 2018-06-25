
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

from esys.escript import *
import esys.finley
from esys.escript.models import StokesProblemCartesian
from esys.finley import finley
from esys.finley import Rectangle
from esys.weipa import saveVTK


#physical properties
rho1 = 1000             #fluid density on bottom
rho2 = 1010             #fluid density on top
eta1 = 100.0            #fluid viscosity on bottom
eta2 = 100.0            #fluid viscosity on top
g=10.0

#solver settings
dt = 0.001
t_step = 0
t_step_end = 2000
TOL = 1.0e-5
max_iter=400
verbose=True
useUzawa=True

#define mesh
l0=0.9142
l1=1.0
n0=200      
n1=200

mesh=Rectangle(l0=l0, l1=l1, order=2, n0=n0, n1=n1)
#get mesh dimensions
numDim = mesh.getDim()
#get element size
h = Lsup(mesh.getSize())

#level set parameters
tolerance = 1.0e-6
reinit_max = 30
reinit_each = 3
alpha = 1
smooth = alpha*h 

#boundary conditions
x = mesh.getX()
#left + bottom + right + top
b_c = whereZero(x[0])*[1.0,0.0] + whereZero(x[1])*[1.0,1.0] + whereZero(x[0]-l0)*[1.0,0.0] + whereZero(x[1]-l1)*[1.0,1.0]

velocity = Vector(0.0, ContinuousFunction(mesh))
pressure = Scalar(0.0, ContinuousFunction(mesh))
Y = Vector(0.0,Function(mesh))

#define initial interface between fluids
xx = mesh.getX()[0]
yy = mesh.getX()[1]
func = Scalar(0.0, ContinuousFunction(mesh))
h_interface = Scalar(0.0, ContinuousFunction(mesh))
h_interface = h_interface + (0.02*cos(math.pi*xx/l0) + 0.2)
func = yy - h_interface
func_new = func.interpolate(ReducedSolution(mesh))

#Stokes Cartesian
solution=StokesProblemCartesian(mesh,debug=True)
solution.setTolerance(TOL)

#level set
levelset = LevelSet(mesh, func_new, reinit_max, reinit_each, tolerance, smooth)    

while t_step <= t_step_end:
  #update density and viscosity
  rho = levelset.update_parameter(rho1, rho2)
  eta = levelset.update_parameter(eta1, eta2)

  #get velocity and pressure of fluid
  Y[1] = -rho*g
  solution.initialize(fixed_u_mask=b_c,eta=eta,f=Y)
  velocity,pressure=solution.solve(velocity,pressure,max_iter=max_iter,verbose=verbose,useUzawa=useUzawa)
  
  #update the interface
  func = levelset.update_phi(velocity, dt, t_step)  

  print("##########################################################")
  print("time step:", t_step, " completed with dt:", dt)
  print("Velocity: min =", inf(velocity), "max =", Lsup(velocity))
  print("##########################################################")
 
  #save interface, velocity and pressure 
  saveVTK("phi2D.%2.4i.vtu"%t_step,interface=func,velocity=velocity,pressure=pressure)
  #courant condition
  dt = 0.4*Lsup(mesh.getSize())/Lsup(velocity)
  t_step += 1
