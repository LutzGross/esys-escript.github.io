
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

############################################################FILE HEADER
# example10a.py
# Model of gravitational Potential for a gravity POLE.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.finley import Rectangle # This imports the rectangle domain function from finley
from esys.weipa import saveVTK # This imports the VTK file saver from weipa
import os, sys #This package is necessary to handle saving our data.
from math import pi, sqrt, sin, cos

from esys.escript.pdetools import Projector, Locator
from cblib import toRegGrid

import matplotlib
matplotlib.use('agg') #It's just here for automated testing

import pylab as pl #Plotting package
import numpy as np

########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    import sys
    print "This example will not run in an MPI world."
    sys.exit(0)

#################################################ESTABLISHING VARIABLES
#Domain related.
mx = 4000*m #meters - model length
my = -4000*m #meters - model width
ndx = 400 # mesh steps in x direction 
ndy = 400 # mesh steps in y direction - one dimension means one element
#PDE related
rho=200.0
rholoc=[mx/2.,my/2.]
G=6.67300*10E-11

R=10
z=50.

################################################ESTABLISHING PARAMETERS
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","example10")
#ensure the dir exists
mkDir(save_path)

#####################################################ANALYTIC SOLUTION
def analytic_gz(x,z,R,drho):
    G=6.67300*10E-11
    return G*2*np.pi*R*R*drho*(z/(x*x+z*z))

sol_angz=[]
sol_anx=[]
for x in range(int(-mx/2),int(mx/2),10):
    sol_angz.append(analytic_gz(x,z,R,rho))
    sol_anx.append(x+mx/2)

##############INVERSION
def gzpot(p, y, x, *args):
    #rho, rhox, rhoy, R = p
    R=10.; rhox=args[0]/2.; rhoy=args[1]/2.
    rho=p[0]
    #Domain related.
    mx = args[0]; my = args[1]; ndx = args[2]; ndy = args[3]

    #PDE related
    G=6.67300*10E-11

    #DOMAIN CONSTRUCTION
    domain = Rectangle(l0=mx,l1=my,n0=ndx, n1=ndy)
    domx=Solution(domain).getX()
    mask=wherePositive(R-length(domx-rholoc))
    rhoe=rho*mask
    kro=kronecker(domain)

    q=whereZero(domx[1]-my)+whereZero(domx[1])+whereZero(domx[0])+whereZero(domx[0]-mx)
    #ESCRIPT PDE CONSTRUCTION
    mypde=LinearPDE(domain)
    mypde.setValue(A=kro,Y=4.*np.pi*G*rhoe,q=q,r=0.0)
    mypde.setSymmetryOn()
    sol=mypde.getSolution()

    g_field=grad(sol) #The graviational accelleration g.
    g_fieldz=g_field*[0,1] #The vertical component of the g field.
    gz=length(g_fieldz) #The magnitude of the vertical component.

    #MODEL SIZE SAMPLING
    sol_escgz=[]
    sol_escx=[]
    for i in range(0,len(x)):
        sol_escgz.append([x[i],my/2.+z])

    sample=[] # array to hold values
    rec=Locator(gz.getFunctionSpace(),sol_escgz) #location to record
    psol=rec.getValue(gz)

    err = np.sum(np.array(y) - np.array(psol))
    print "Lsup= ",Lsup(np.array(psol)-np.array(sol_angz))/Lsup(np.array(psol))
    return err

#Initial Guess
#guess=[400,mx/4,my/4,50]
guess=15.

from scipy.optimize import leastsq
#plsq = leastsq(gzpot, guess, args=(sol_angz, sol_anx, mx, my, ndx, ndy),maxfev=10)

objf=[]
for R in range(180,220):
    objf.append(gzpot([R],sol_angz,sol_anx, mx, my, ndx, ndy))

pl.plot(objf)
pl.show()

