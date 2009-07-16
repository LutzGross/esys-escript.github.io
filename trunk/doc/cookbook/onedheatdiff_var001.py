
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



##ESTABLISHING VARIABLES
#PDE related
mx = 500 # model lenght
my = 100 # model width
ndx = 500 # steps in x direction 
ndy = 1 # steps in y direction

q=0 #our heat source temperature is now zero
Tref=2273 # Kelvin #the starting temperature of our intrusion
rho = 2750 #kg/m^{3} density
cp = 790 #j/kg specific heat
rhocp = rho*cp	#DENSITY * SPECIFIC HEAT
eta=0.  # RADIATION CONDITION
kappa=2.2 # Watts/(meter*Kelvin) DIFFUSION CONSTANT, HEAT PERMEABILITY
#Script/Iteration Related
t=0. #our start time, usually zero
tday=10*365. #the time we want to end the simulation in days
tend=tday*24*60*60
outputs = 400 # number of time steps required.
h=(tend-t)/outputs #size of time step

print "Expected Number of Output Files is: ", outputs
print "Step size is: ", h/(24.*60*60), "days"


i=0 #loop counter 
save_path = "data/onedheatdiff_var001" #the folder to put our outputs in, leave blank "" for script path - note this folder path must exist to work

#... generate domain ...
model = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
# extract finite points
x=model.getX()
#... open PDE ...
mypde=LinearPDE(model)
mypde.setSymmetryOn()
mypde.setValue(A=kappa*kronecker(model),D=rhocp/h,d=eta,y=eta*Tref)

# ... set initial temperature ....
bound = x[0]-mx/(ndx/250.)
T= 0*Tref*whereNegative(bound)+Tref*wherePositive(bound)

saveVTK(os.path.join(save_path,"data%03d.vtu") %i,sol=T)

# ... start iteration:
while t<=tend:
      i+=1
      t+=h
      mypde.setValue(Y=rhocp/h*T)
      T=mypde.getSolution()
      saveVTK(os.path.join(save_path,"data%03d.vtu") %i,sol=T)
      
# iteration var 2
#Tl = 0
#Tr = Tref
#while Tl < Tr*0.8:
      #mypde.setValue(Y=rhocp/h*T)
      #T=mypde.getSolution()
      #i+=1
      #x=rod.getX()
      #Tl= x[0]
      #Tr= x[

#print "Finish temp balance in:", i*h/(24.*60*60), " days"
