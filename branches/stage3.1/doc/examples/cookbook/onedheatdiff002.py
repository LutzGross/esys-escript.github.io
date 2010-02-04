
########################################################
#
# Copyright (c) 2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009 by University of Queensland
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
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.finley import Rectangle # This imports the rectangle domain function from finley
#For interactive use, you can comment out the next two lines
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
import pylab as pl #Plotting package.
import numpy as np #Array package.
import os, sys #This package is necessary to handle saving our data.

# .. MPI WORLD CHECK
if getMPISizeWorld() > 1:
    import sys
    print "This example will not run in an MPI world."
    sys.exit(0)

##ESTABLISHING VARIABLES
#Domain related.
mx = 1*m #meters - model length
my = .1*m #meters - model width
ndx = 100 # mesh steps in x direction 
ndy = 1 # mesh steps in y direction - one dimension means one element
#PDE related
rho = 7874. *kg/m**3 #kg/m^{3} density of iron
cp = 449.*J/(kg*K) # J/Kg.K thermal capacity
rhocp = rho*cp 
kappa = 80.*W/m/K   # watts/m.Kthermal conductivity
qH=0 * J/(sec*m**3) # J/(sec.m^{3}) no heat source
Tref = 20 * Celsius  # base temperature of the rod
T0 = 100 * Celsius # temperature at heating element

t=0 * day  # our start time, usually zero
tend= 0.5 *day  # - time to end simulation
outputs = 200 # number of time steps required.
h=(tend-t)/outputs #size of time step
#user warning statement
print "Expected Number of time outputs is: ", (tend-t)/h
i=0 #loop counter
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","onedheatdiff002")
#ensure the dir exists
mkDir(save_path, os.path.join(save_path,"tempT"))

#... generate domain ...
rod = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
x=Solution(rod).getX()
#... open PDE and set coefficients ...
mypde=LinearPDE(rod)
A=zeros((2,2))
A[0,0]=kappa
q=whereZero(x[0])
mypde.setValue(A=A, D=rhocp/h, q=q, r=T0)
# ... set initial temperature ....
T= T0*whereZero(x[0])+Tref*(1-whereZero(x[0]))

# ... open a collector for the time marks and corresponding total energy
t_list=[]
E_list=[]
# ... convert solution points for plotting
plx = x.toListOfTuples() 
plx = np.array(plx) #convert to tuple to numpy array
plx = plx[:,0] #extract x locations
# ... start iteration:
while t<tend:
      i+=1
      t+=h
      mypde.setValue(Y=qH+rhocp/h*T)
      T=mypde.getSolution()
      totE=integrate(rhocp*T)
      print "time step %s at t=%e minutes completed. total energy = %e."%(i,t/minute,totE)
      t_list.append(t)
      E_list.append(totE)

      #establish figure 1 for temperature vs x plots
      tempT = T.toListOfTuples()
      pl.figure(1) #current figure
      pl.plot(plx,tempT) #plot solution
      # add title
      pl.axis([0,mx,Tref*.9,T0*1.1])
      pl.title("Temperature across rod at time %e hours"%(t/hour))
      #save figure to file
      pl.savefig(os.path.join(save_path,"tempT", "rodpyplot%03d.png"%i))
      pl.clf() #clear figure
# plot the total energy over time:
pl.figure(2)
pl.plot(t_list,E_list)
pl.title("Total Energy")
# pl.axis([0,max(t_list),0,max(E_list)*1.1])
pl.savefig(os.path.join(save_path,"totE.png"))
pl.clf()

# compile the *.png files to create a*.avi video that show T change
# with time. This opperation uses linux mencoder. For other operating 
# systems it may be possible to use your favourite video compiler to
# convert image files to videos. To enable this step uncomment the
# following lines.

#os.system("mencoder mf://"+save_path+"/tempT"+"/*.png -mf type=png:\
#w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
#onedheatdiff001tempT.avi")
