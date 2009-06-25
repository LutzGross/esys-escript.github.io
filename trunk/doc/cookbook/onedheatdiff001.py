
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
from esys.escript.unitsSI import * # A useful unit handling package which will make sure all our units match up in the equations.
import os #This package is necessary to handle saving our data.
#plotting tools
#import matplotlib as ptool


##ESTABLISHING VARIABLES
#Domain related.
mx = 1*m #meters - model lenght
my = .1*m #meters - model width
ndx = 100 # steps in x direction 
ndy = 1 # steps in y direction

#PDE related
q=200. * Celsius #Kelvin - our heat source temperature
Tref = 0. * Celsius # Kelvin - starting temp of iron bar
rho = 7874. *kg/m**3 #kg/m^{3} density of iron
cp = 449.*J/(kg*K) #j/Kg.K thermal capacity
rhocp = rho*cp
kappa = 80.*W/m/K #watts/m.Kthermal conductivity
#Script/Iteration Related
t=0 #our start time, usually zero
tend=5.*minute #seconds - time to end simulation
outputs = 200 # number of time steps required.
h=(tend-t)/outputs #size of time step
print "Expected Number of time outputs is: ", (tend-t)/h
i=0 #loop counter
#the folder to put our outputs in, leave blank "" for script path 
save_path="data/onedheatdiff001"
#note this folder path must exist to work 


#... generate domain ...
rod = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
# extract finite points
x=rod.getX()
#... open PDE ...
mypde=LinearPDE(rod)
mypde.setSymmetryOn()
mypde.setValue(A=kappa*kronecker(rod),D=rhocp/h)

# ... set heat source: ....
qH=q*whereZero(x[0])
# ... set initial temperature ....
T=Tref

# ... start iteration:
while t<=tend:
      i+=1
      t+=h
      mypde.setValue(Y=qH+rhocp/h*T)
      T=mypde.getSolution()
      #ptool.plot(T)
      saveVTK(os.path.join(save_path,"data%03d.xml") %i,sol=T)
      

