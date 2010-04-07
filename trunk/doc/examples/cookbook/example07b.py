
########################################################
#
# Copyright (c) 2009-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# Antony Hallam
# Acoustic Wave Equation Simulation using acceleration solution
# and lumping.

# Importing all the necessary modules required.
from esys.escript import *
from esys.finley import Rectangle
import sys
import os
# smoothing operator 
from esys.escript.pdetools import Projector
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from esys.escript.linearPDEs import LinearPDE

# Establish a save path.
savepath = "data/example07b"
mkDir(savepath)

#Geometric and material property related variables.
mx = 1000. # model lenght
my = 1000. # model width
ndx = 400 # steps in x direction 
ndy = 400 # steps in y direction

xstep=mx/ndx
ystep=my/ndy

c=380.0
csq=c*c
# Time related variables.
tend=1.5    #end time
#calculating )the timestep
h=0.001
#Check to make sure number of time steps is not too large.
print "Time step size= ",h, "Expected number of outputs= ",tend/h

#uncomment the following lines to give the user a chance to stop
#proceeder = raw_input("Is this ok?(y/n)")
#Exit if user thinks too many outputs.
#if proceeder == "n":
#   sys.exit()

U0=0.01 # amplitude of point source
#  spherical source at middle of bottom face

xc=[500,500]

mydomain=Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
#wavesolver2d(mydomain,h,tend,lam,mu,rho,U0,xc,savepath,output="mpl")
x=mydomain.getX()

# ... open new PDE ...
mypde=LinearPDE(mydomain)
print mypde.isUsingLumping()
print mypde.getSolverOptions()
mypde.getSolverOptions().setSolverMethod(mypde.getSolverOptions().LUMPING)
mypde.setSymmetryOn()
mypde.setValue(D=1.)

# define small radius around point xc
src_radius = 30
print "src_radius = ",src_radius

# ... set initial values ....
n=0
# for first two time steps
u=U0*(cos(length(x-xc)*3.1415/src_radius)+1)*whereNegative(length(x-xc)-src_radius)
u_m1=u
t=0

#plot source shape
uT=np.array(u.toListOfTuples())
uT=np.reshape(uT,(ndx+1,ndy+1))
source_line=uT[ndx/2,:]
pl.plot(source_line)
pl.plot(source_line,'ro')
pl.axis([70,130,0,0.05])
pl.savefig(os.path.join(savepath,"source_line.png"))

while t<tend:
    # get current pressure
    g=grad(u)
    pres=csq*g
    # set values
    mypde.setValue(X=-pres)
    # get new acceleration
    accel = mypde.getSolution()
    # shift displacements
    u_p1=(2.*u-u_m1)+h*h*accel
    u_m1=u
    u=u_p1
	# increment loop values
    t+=h
    n+=1
    print n,"-th time step t ",t
    # ... save current acceleration in units of gravity and displacements 
    saveVTK(os.path.join(savepath,"ex07b.%i.vtu"%n),output1 = length(u),tensor=pres)

#~ u_pc_data.close()
#~ os.system("mencoder mf://"+savepath+"/*.png -mf type=png:\
#~ w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
#~ wsmpl.avi")

#mencoder mf://*.png -mf type=png:\w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o wsmpl.avi
