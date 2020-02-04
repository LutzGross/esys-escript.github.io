__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

"""
Simple direct current resistivity code that uses class dcResistivityModel in dcModels.py
The domain is made using finley Brick creating a structured grid.  
Inputs are 
    gridline spacing,
    number of elements in the x, y, and z directions,
    height of the data above ground,
    background magnetic field,
    and assumed susceptibility.

class dcModel(domain) sets up the PDE and domain

This example is for a model with conductivity = 1 apart from a small sphere centre c and radius R that has the assumed conductivity.

Sources are on the top of the domain, z=zmax.

The output from this code is a silo file containing conductivity and primary and secondary electric potential.   
"""

# Import required modules
from esys.escript import *
from esys.finley import Brick
from esys.weipa import saveSilo
from esys.downunder.apps import DCResistivityModel
import numpy as np 

# Set Parameters
dx = 6         # grid line spacing in [m]
NEx = 200      # number of nodes in the x direction
NEy = 200      # number of nodes in the y direction
NEz = 200      # number of nodes in the z direction
H0 = 600       # height [m] of transect above bottom of domain (will be locked to grid)
sig_p = 1.     # primary conductivity
sig_2 = 100.   # conductivity in ball
useAnalytic = True

POSx = 100
POSy = 70

NEGx = 100
NEGy = 110

POS_node = ( POSx*dx, POSy*dx, NEz*dx)
NEG_node = ( NEGx*dx, NEGy*dx, NEz*dx)

Lx = dx*NEx
Ly = dx*NEy
Lz = dx*NEz
print("domain dimensions = [%dm x %dm x %dm]"%(Lx,Ly,Lz))
print("grid = [%d x %d x %d]"%(NEx, NEy, NEz))

# Create domain with Dirac points
domain=Brick(n0=NEx, n1=NEy, n2=NEz, l0=Lx, l1=Ly, l2=Lz,
              diracPoints= [POS_node, NEG_node], diracTags = ['e0','e1'])

# Define model using class dcresistivityModel3D
# Assumes zero Dirichlet BC at bottom left front corner


# Define conductivity 
Z0=100           # vertical position of circle below transect [m]
c=[Lx/2.,Ly/2., H0-Z0] # circle center
R=50.            # radius

x=domain.getX()
d=length(x-c)

sigma = Scalar(sig_p,ContinuousFunction(domain))
model=DCResistivityModel(domain,sigma0=sigma,useFastSolver = True)
print("ERT forward model created.")
print("background conductivity = %s"%(sigma))

print("positive charge %s at %s"%(1.0, POS_node))
print("negative charge %s at %s"%(1.0, NEG_node))

if useAnalytic:
    print("half space solution as primary potential used.")
    model.setPrimaryPotentialForHalfSpace(sources= [POS_node, NEG_node], 
                                      charges=[1.0, -1.0] ) 
    print("primary potential set.")
else:
    src=Scalar(0., DiracDeltaFunctions(domain))
    src.setTaggedValue('e0', 1.0)
    src.setTaggedValue('e1', -1.0)
    print("source = %s"%src)
    model.setPrimaryPotential(source=src) 
    print("primary potential set.")    

sphereCond=sigma+(sig_2-sig_p)*whereNegative(d-R)    # 0 for d>R and 1 for d<R

model.setConductivity(sphereCond)

saveSilo("dcAppExample",sigma=model.getConductivity(), u2=model.getSecondaryPotential(),up=model.getPrimaryPotential())

