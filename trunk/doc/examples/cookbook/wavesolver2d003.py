##############################################################################
#
# Copyright (c) 2009-2018 by The University of Queensland
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
from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# You can shorten the execution time by reducing variable tend from 60 to 0.5
# Antony Hallam
# Acoustic Wave Equation Simulation

# Importing all the necessary modules required.
import matplotlib
matplotlib.use('agg') #It's just here for automated testing
import os
import sys
import numpy as np
import pylab as pl
import matplotlib.cm as cm

from esys.escript import *
# smoothing operator 
from esys.escript.pdetools import Projector
from esys.finley import Rectangle
from esys.weipa import saveVTK
from cblib1 import wavesolver2d

# Establish a save path.
savepath = "data/wavesolver2d009mpltestABCnolump0_0006"
mkDir(savepath)


#Geometric and material property related variables.
mx = 1000. # model lenght
my = 1000. # model width
ndx = 200 # steps in x direction 
ndy = 200 # steps in y direction

xstep=mx/ndx
ystep=my/ndy

lam=3.462e9 #lames constant
mu=3.462e9  #bulk modulus
rho=1154.   #density
# Time related variables.
tend=0.5    #end time
#calculating )the timestep
h=(1./5.)*sqrt(rho/(lam+2*mu))*(mx/ndx)
#Check to make sure number of time steps is not too large.
print("Time step size= ",h, "Expected number of outputs= ",tend/h)

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




domain=mydomain
output="mpl"





from esys.escript.linearPDEs import LinearPDE, SolverOptions
x=domain.getX()

## boundary conditions

bleft=xstep*50.
bright=mx-(xstep*50.)
bbot=my-(ystep*50.)
btop=ystep*50.

left=x[0]-bleft
right=x[0]-bright
bottom=x[1]-bbot
top=x[1]-btop

decay=0.0006
fleft=exp(-1.0*(decay*(bleft-x[0]))**2)
fright=exp(-1.0*(decay*(x[0]-bright))**2)
fbottom=exp(-1.0*(decay*(x[1]-bbot))**2)
ftop=exp(-1.0*(decay*(btop-x[1]))**2)

abcleft=fleft*whereNegative(left)
abcright=fright*wherePositive(right)
abcbottom=fbottom*wherePositive(bottom)
abctop=ftop*whereNegative(top)

abcleft=abcleft+whereZero(abcleft)
abcright=abcright+whereZero(abcright)
abcbottom=abcbottom+whereZero(abcbottom)
abctop=abctop+whereZero(abctop)

abc=abcleft*abcright*abcbottom*abctop

#~ fleftT=fleft.toListOfTuples()
#~ fleftT=np.reshape(fleftT,(ndx+1,ndy+1))
#~ pl.imshow(fleftT)
#~ pl.colorbar()
#~ pl.savefig("fleftT.png")
#~ 
#~ frightT=fright.toListOfTuples()
#~ frightT=np.reshape(frightT,(ndx+1,ndy+1))
#~ pl.clf()
#~ pl.imshow(frightT)
#~ pl.colorbar()
#~ pl.savefig("frightT.png")
#~ 
#~ fbottomT=fbottom.toListOfTuples()
#~ fbottomT=np.reshape(fbottomT,(ndx+1,ndy+1))
#~ pl.clf()
#~ pl.imshow(fbottomT)
#~ pl.colorbar()
#~ pl.savefig("fbottomT.png")
#~ 
#~ #tester=fright*wherePositive(right)
#~ tester=fleft*whereNegative(left)
#~ tester=tester.toListOfTuples()
#~ tester=np.reshape(tester,(ndx+1,ndy+1))
#~ pl.clf()
#~ pl.imshow(tester)
#~ pl.colorbar()
#~ pl.savefig("tester1.png")
#~ 
#~ tester=fright*wherePositive(right)
#~ tester=tester.toListOfTuples()
#~ tester=np.reshape(tester,(ndx+1,ndy+1))
#~ pl.clf()
#~ pl.imshow(tester)
#~ pl.colorbar()
#~ pl.savefig("tester2.png")
#~ 
#~ tester=fbottom*wherePositive(bottom)
#~ tester=tester.toListOfTuples()
#~ tester=np.reshape(tester,(ndx+1,ndy+1))
#~ pl.clf()
#~ pl.imshow(tester)
#~ pl.colorbar()
#~ pl.savefig("tester3.png")


# ... open new PDE ...
mypde=LinearPDE(domain)
print(mypde.isUsingLumping())
print(mypde.getSolverOptions())
#mypde.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
mypde.setSymmetryOn()
kmat = kronecker(domain)
mypde.setValue(D=kmat*rho)

# define small radius around point xc
# Lsup(x) returns the maximum value of the argument x
src_radius = 50#2*Lsup(domain.getSize())
print("src_radius = ",src_radius)

dunit=numpy.array([0.,1.]) # defines direction of point source
#~ dunit=(x-xc)
#~ absrc=length(dunit)
#~ dunit=dunit/maximum(absrc,1e-10)

# ... set initial values ....
n=0
# initial value of displacement at point source is constant (U0=0.01)
# for first two time steps
u=U0*(cos(length(x-xc)*3.1415/src_radius)+1)*whereNegative(length(x-xc)-src_radius)*dunit
#u=whereNegative(length(x-xc)-src_radius)*dunit

maxi=0.02

print(u)
u_m1=u
t=0

#~ u_pot = cbphones(domain,u,[[0,500],[250,500],[400,500]],2)
#~ u_pc_x1 = u_pot[0,0]
#~ u_pc_y1 = u_pot[0,1]
#~ u_pc_x2 = u_pot[1,0]
#~ u_pc_y2 = u_pot[1,1]
#~ u_pc_x3 = u_pot[2,0]
#~ u_pc_y3 = u_pot[2,1]
#~ 
#~ # open file to save displacement at point source
#~ u_pc_data=open(os.path.join(savepath,'U_pc.out'),'w')
#~ u_pc_data.write("%f %f %f %f %f %f %f\n"%(t,u_pc_x1,u_pc_y1,u_pc_x2,u_pc_y2,u_pc_x3,u_pc_y3))

while t<tend:
    # ... get current stress ....
#    t=1.
    ##OLD WAY
    g=grad(u)
    stress=lam*trace(g)*kmat+mu*(g+transpose(g))
    ### ... get new acceleration ....
    #mypde.setValue(X=-stress)          
    #a=mypde.getSolution()
    ### ... get new displacement ...
    #u_p1=2*u-u_m1+h*h*a
    ###NEW WAY
    mypde.setValue(X=-stress*(h*h),Y=(rho*2*u-rho*u_m1))
    u_p1 = mypde.getSolution()
    # ... shift displacements ....
    u_m1=u
    u=u_p1*abc
    #stress = 
    t+=h
    n+=1
    print(n,"-th time step t ",t)
    #~ u_pot = cbphones(domain,u,[[300.,200.],[500.,200.],[750.,200.]],2)
    #~ 
    #~ #     print "u at point charge=",u_pc
    #~ u_pc_x1 = u_pot[0,0]
    #~ u_pc_y1 = u_pot[0,1]
    #~ u_pc_x2 = u_pot[1,0]
    #~ u_pc_y2 = u_pot[1,1]
    #~ u_pc_x3 = u_pot[2,0]
    #~ u_pc_y3 = u_pot[2,1]

    # save displacements at point source to file for t > 0
    #~ u_pc_data.write("%f %f %f %f %f %f %f\n"%(t,u_pc_x1,u_pc_y1,u_pc_x2,u_pc_y2,u_pc_x3,u_pc_y3))

    # ... save current acceleration in units of gravity and displacements 
    #saveVTK(os.path.join(savepath,"usoln.%i.vtu"%n),acceleration=length(a)/9.81,
    #displacement = length(u), tensor = stress, Ux = u[0] )
    if output == "vtk":
        saveVTK(os.path.join(savepath,"tonysol.%i.vtu"%n),output1 = length(u),tensor=stress)
    if output == "mpl":
        uT=np.array(u.toListOfTuples())
        uT=np.reshape(uT,(ndx+1,ndy+1,2))
        uTz=uT[:,:,1]+uT[:,:,0]
        uTz=np.transpose(uTz)
        pl.clf()
        # plot wave
        uTz[0,0]=maxi
        uTz[0,1]=-maxi
        CS = pl.imshow(uTz,cmap=cm.spectral)
        pl.colorbar()
        # labels and formatting
        pl.title("Wave Equation Cookbook Example ABC.")
        pl.xlabel("Horizontal Displacement (m)")
        pl.ylabel("Depth (m)")
        if getMPIRankWorld() == 0: #check for MPI processing
            pl.savefig(os.path.join(savepath,"ws04mpl%05d.png"%n))

#~ u_pc_data.close()
#~ os.system("mencoder mf://"+savepath+"/*.png -mf type=png:\
#~ w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o \
#~ wsmpl.avi")

#mencoder mf://*.png -mf type=png:\w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o wsmpl.avi
