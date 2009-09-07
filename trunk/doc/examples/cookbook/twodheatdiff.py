
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.finley import Rectangle # This imports the rectangle domain function from finley
import os #This package is necessary to handle saving our data.

from cblib import needdirs


##ESTABLISHING VARIABLES
#PDE related
mx = 600 # model lenght
my = 600 # model width
ndx = 100 # steps in x direction 
ndy = 100 # steps in y direction
r = 200 # radius of intrusion
ic = [300, 0] #centre of intrusion

q=0 #our heat source temperature is now zero
Ti=2273 # Kelvin #the starting temperature of our iron bar
rhoi = 2750 #kg/m^{3} density
cpi = 790 #j/kg specific heat
rhocpi = rhoi*cpi	#DENSITY * SPECIFIC HEAT
eta=0.  # RADIATION CONDITION
kappai=2.2 # Watts/(meter*Kelvin) DIFFUSION CONSTANT, HEAT PERMEABILITY

Tc = 200
rhoc = 2200
cpc = 400
rhocpc = rhoc*cpc
kappac = 0.1


#Script/Iteration Related
t=0. #our start time, usually zero
tday=100*365. #the time we want to end the simulation in days
tend=tday*24*60*60
outputs = 200 # number of time steps required.
h=(tend-t)/outputs #size of time step

print "Expected Number of Output Files is: ", outputs
print "Step size is: ", h/(24.*60*60), "days"


i=0 #loop counter 
save_path = "data/twodheatdiff" #the folder to put our outputs in, leave blank "" for script path - note this folder path must exist to work
needdirs([save_path])

#... generate domain ...
model = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
# extract finite points
x=model.getX()

#... open PDE ...
mypde=LinearPDE(model)
mypde.setSymmetryOn()

bound = length(x-ic)-r #where the boundary will be located

A = (kappai)*whereNegative(bound)+(kappac)*wherePositive(bound)
D = (rhocpi/h)*whereNegative(bound)+(rhocpc/h)*wherePositive(bound)

mypde.setValue(A=A*kronecker(model),D=D,d=eta,y=eta*Tc)

# ... set initial temperature ....

T= Ti*whereNegative(bound)+Tc*wherePositive(bound) #defining the initial temperatures.
saveVTK(os.path.join(save_path,"dataedge.vtu"), sol=bound)
saveVTK(os.path.join(save_path,"data%03d.vtu") %i,sol=T)

#... start iteration:
while t<=tend:
      i+=1
      t+=h
      Y = T*D
      mypde.setValue(Y=Y)
      T=mypde.getSolution()
      saveVTK(os.path.join(save_path,"data%03d.vtu") %i,sol=T)
