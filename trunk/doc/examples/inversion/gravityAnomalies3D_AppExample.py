__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

"""
Simple gravity Anomaly code that uses class GravityModel in gravityModels.py
The domain is made using finley Brick creating a structured grid.  
Inputs are 
    gridline spacing,
    number of elements in the x, y, and z directions,
    height of the data above ground,
    background magnetic field,
    and assumed susceptibility.

class GravityModel(domain) sets up the PDE and domain

This example is for a model with zero density apart from a small sphere centre c and radius R that has the assumed density.

The output from this code is a silo file containing density, and vertical gravity, gz.   
"""

# Import required modules
from esys.escript import *
from esys.finley import Brick
from esys.weipa import saveSilo
from esys.downunder.apps import GravityModel
import numpy as np 

# Set Parameters
dx = 6         # grid line spacing in [m]
NEx = 200      # number of nodes in the x direction
NEy = 200      # number of nodes in the y direction
NEz = 200      # number of nodes in the z direction
H0 = 600       # height [m] of transect above bottom of domain (will be locked to grid)
dens = 1000     # assumed density kg/m^3


Lx = dx*NEx
Ly = dx*NEy
Lz = dx*NEz
print("domain dimensions = [%dm x %dm x %dm]"%(Lx,Ly,Lz))
print("grid = [%d x %d x %d]"%(NEx, NEy, NEz))

# Create domain
domain=Brick(n0=NEx, n1=NEy, n2=NEz, l0=Lx, l1=Ly, l2=Lz)

# Define model using class GravityModel
# Assumes zero Dirichlet BC at top surface
model=GravityModel(domain)

# Define density 
Z0=100           # vertical position of circle below transect [m]
c=[Lx/2.,Ly/2., H0-Z0] # circle center
R=50.            # radius

x=domain.getX()
d=length(x-c)
sphereDens=dens*whereNegative(d-R)    # 0 for d>R and 1 for d<R
model.setDensity(sphereDens)

# Solve for gravity field anomaly
gp=model.getGravityPotential()
print(gp)

g_vect = model.getGravityVector()
gz = model.getzGravity()
rho=model.getDensity()
print(rho)

saveSilo("grav", gz=gz, rho=rho ,allg=g_vect)
