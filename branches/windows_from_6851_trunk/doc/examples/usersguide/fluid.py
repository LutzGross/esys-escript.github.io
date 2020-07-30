##############################################################################
#
# Copyright (c) 2008-2018 by The University of Queensland
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

########      August 2008      ########
##########    Leon Graham    ########## 
## Newtonian fluid using StokesProblemCartesian class##
from __future__ import division, print_function

from esys.escript import *
try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError:
    print("Finley module required but not available")
    HAVE_FINLEY = False
from esys.escript.linearPDEs import LinearPDE
from esys.escript.models import StokesProblemCartesian
from esys.weipa import saveVTK
if HAVE_FINLEY:
    #physical constants
    eta=1.0
    rho=100.0
    g=10.0 

    #solver settings
    tolerance=1.0e-4
    max_iter=200
    t_end=50
    t=0.0
    time=0
    verbose='TRUE'
    useUzawa='TRUE'

    #define mesh 
    H=2.0
    L=1.0
    W=1.0
    mesh = Rectangle(l0=L, l1=H, order=-1, n0=20, n1=20, useElementsOnFace=0) # use linear macro elements for pressure
    coordinates = mesh.getX()

    #gravitational force
    Y=Vector(0.0, Function(mesh))
    Y[1]=-rho*g

    #element spacing
    h=Lsup(mesh.getSize())

    #boundary conditions for slip at base
    boundary_cond=whereZero(coordinates[1])*[0.0,1.0]+whereZero(coordinates[0])*[1.0,0.0]

    #velocity and pressure vectors
    velocity=Vector(0.0, Solution(mesh))
    pressure=Scalar(0.0, ReducedSolution(mesh))

    #Stokes Cartesian
    solution=StokesProblemCartesian(mesh)
    solution.setTolerance(tolerance)

    while t <= t_end:

      print(" ----- Time step = %s -----"%( t ))
      print("Time = %s seconds"%( time ))  
     
      solution.initialize(fixed_u_mask=boundary_cond,eta=eta,f=Y)
      velocity,pressure=solution.solve(velocity,pressure,max_iter=max_iter,verbose=verbose,usePCG=True)
      
      print("Max velocity =", Lsup(velocity), "m/s")
      
      #Courant condition
      dt=0.4*h/(Lsup(velocity))
      print("dt", dt) 
      
      #displace the mesh
      displacement = velocity * dt
      coordinates = mesh.getX()
      newx=interpolate(coordinates + displacement, ContinuousFunction(mesh))
      mesh.setX(newx)  
      
      time += dt
      
      vel_mag = length(velocity)

      #save velocity and pressure output
      saveVTK("vel.%2.2i.vtu"%(t),vel=vel_mag,vec=velocity,pressure=pressure)
      t = t+1.0
