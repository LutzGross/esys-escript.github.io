__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

"""
Simple magnetic Anomaly code that uses class MagneticModel3D in magneticModels.py
The domain is made using finley Brick creating a structured grid.  
Inputs are 
    gridline spacing,
    number of elements in the x, y, and z directions,
    height of the data above ground,
    background magnetic field,
    and assumed susceptibility.

class MagneticModel3D(domain) sets up the PDE and domain

This example is for a model with zero magnetic susceptibility apart from a small sphere centre c and radius R that has the assumed susceptibility.

The output from this code is a silo file containing susceptibility, k, and magnetic field, ba.   
"""

# Import required modules
from esys.escript import *
from esys.finley import Brick
from esys.weipa import saveSilo
from esys.downunder.apps import MagneticModel3D
import numpy as np 

# Set Parameters
dx = 6         # grid line spacing in [m]
NEx = 200      # number of nodes in the x direction
NEy = 200      # number of nodes in the y direction
NEz = 200      # number of nodes in the z direction
H0 = 600       # height [m] of transect above bottom of domain (will be locked to grid)
b_hx = 45000.0 # background magnetic field in nT x direction
b_hz = 0.0     # background magnetic field in nT y direction
b_hy = 0.0     # background magnetic field in nT z direction
ksi = 0.015    # assumed susceptibility


Lx = dx*NEx
Ly = dx*NEy
Lz = dx*NEz
print("domain dimensions = [%dm x %dm x %dm]"%(Lx,Ly,Lz))
print("grid = [%d x %d x %d]"%(NEx, NEy, NEz))

# Create domain
domain=Brick(n0=NEx, n1=NEy, n2=NEz, l0=Lx, l1=Ly, l2=Lz)

# Define model using class MagneticModel3D
# Assumes zero Dirichlet BC at bottom left front corner
model=MagneticModel3D(domain)

# Define susceptibility 
Z0=100           # vertical position of circle below transect [m]
c=[Lx/2.,Ly/2., H0-Z0] # circle center
R=20.            # radius

x=domain.getX()
d=length(x-c)
kC=ksi*whereNegative(d-R)    # 0 for d>R and 1 for d<R
model.setSusceptibility(kC)

# Set background magnetic field [Bx, By]
B_h=[-b_hx, -b_hy, -b_hz]
model.setBackgroundMagneticField(B_h)

# Solve for magnetic field anomaly
B_a=model.getMagneticFieldAnomaly()
print(B_a)
b_a=length(B_a+B_h)-length(np.array(B_h))

#B_h=[0, -b_h]
#model.setBackgroundMagneticField(B_h)

#B_a=model.getMagneticFieldAnomaly()

#b_a=length(B_a+B_h)-b_h

saveSilo("mag3D", ba=b_a, k=model.getSusceptibility())
