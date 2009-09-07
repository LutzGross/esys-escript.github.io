
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
q=50.e6 #our heat source temperature
Tref=0. #the starting temperature of our iron bar
rho=2.6e6
eta=0#75.
kappa=240.
#Script/Iteration Related
t=0 #our start time, usually zero
tend=5.#the time we want to end the simulation
h=0.05 #size of time step

print "Expected Number of Output Files is: ", (tend-t)/h

i=0 #loop counter 
save_path = "data/onedheatdiff" #the folder to put our outputs in, leave blank "" for script path - note this folder path must exist to work
needdirs([save_path])

#... generate domain ...
rod = Rectangle(l0=0.05,l1=.01,n0=500, n1=1)
# extract finite points
x=rod.getX()
#... open PDE ...
mypde=LinearPDE(rod)
mypde.setSymmetryOn()
mypde.setValue(A=kappa*kronecker(rod),D=rho/h,d=eta,y=eta*Tref)
# ... set heat source: ....

qH=q*whereZero(x[0])
# ... set initial temperature ....
T=Tref

# ... start iteration:
while t<=tend:
      i+=1
      t+=h
      mypde.setValue(Y=qH+rho/h*T)
      T=mypde.getSolution()
      print T
      saveVTK(os.path.join(save_path,"data%03d.vtu") %i,sol=T)
      
#~ command = ('mencoder',
           #~ 'mf://*.png',
           #~ '-mf',
           #~ 'type=png:w=800:h=600:fps=25',
           #~ '-ovc',
           #~ 'lavc',
           #~ '-lavcopts',
           #~ 'vcodec=mpeg4',
           #~ '-oac',
           #~ 'copy',
           #~ '-o',
           #~ 'output.avi')

