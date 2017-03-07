##############################################################################
#
# Copyright (c) 2009-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2009-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

############################################################FILE HEADER
# example10d.py
# Model of gravitational Potential for a gravity POLE.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
import matplotlib
matplotlib.use('agg') #Its just here for automated testing

from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.weipa import saveVTK # This imports the VTK file saver from weipa
import os, sys #This package is necessary to handle saving our data.
from math import pi, sqrt, sin, cos

from esys.escript.pdetools import Projector, Locator
from esys.finley import ReadGmsh

import pylab as pl #Plotting package
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

try:
    # This imports the rectangle domain function 
    from esys.finley import MakeDomain
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
    import sys
    print("This example will not run in an MPI world.")
    sys.exit(0)

if HAVE_FINLEY:
    #################################################ESTABLISHING VARIABLES
    #Domain related.
    mx = 10000*m #meters - model length
    my = 10000*m #meters - model width
    #PDE related
    rho=10.0
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
    for x in range(int(-mx/20),int(mx/20),10):
        sol_angz.append(analytic_gz(x,z,R,rho))
        sol_anx.append(x+mx/2)

    ##############INVERSION
    def gzpot(p, y, x, *args):
        #rho, rhox, rhoy, R = p
        rhox=args[0]/2.; rhoy=args[1]/2.
        rho, R, z =p
        #Domain related.
        mx = args[0]; my = args[1];
        #PDE related
        G=6.67300*10E-11

        #DOMAIN CONSTRUCTION
        domain=ReadGmsh('data/example10m/example10m.msh',2)
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
            sol_escgz.append([x[i],rhoy+z])

        sample=[] # array to hold values
        rec=Locator(gz.getFunctionSpace(),sol_escgz) #location to record
        psol=rec.getValue(gz)

        err = np.sum((np.array(y) - np.array(psol))**2.)
        print("Lsup= ",Lsup(np.array(psol)-np.array(sol_angz))/Lsup(np.array(psol)))
        return err

    #Initial Guess
    #guess=[400,mx/4,my/4,50]
    guess=[15.,20.]

    from scipy.optimize import leastsq
    #plsq = leastsq(gzpot, guess, args=(sol_angz, sol_anx, mx, my, ndx, ndy),maxfev=20)
    #print plsq

    objf=[]
    x=np.arange(1,20)
    y=np.arange(1,20)
    z=np.arange(40,60)
    fig=pl.figure(figsize=(5,5))

    for p in x:
        objf.append(gzpot([p,10.,50.],sol_angz,sol_anx, mx, my))
    sp=fig.add_subplot(311)
    sp.plot(x,objf)
    sp.set_title("Variable RHO")
    objf=[]
    for R in y:
        objf.append(gzpot([10.,R,50.],sol_angz,sol_anx, mx, my))
    sp=fig.add_subplot(312)
    sp.plot(y,objf)
    sp.set_title("Variable Radius")

    objf=[]
    for D in z:
        objf.append(gzpot([10.,10.,D],sol_angz,sol_anx, mx, my))     
    sp=fig.add_subplot(313)
    sp.plot(z,objf)
    sp.set_title("Variable Depth")

    fig.savefig("ex10e_objf.pdf",dpi=150)

    #ob=np.array(objf)
    #X,Y=pl.meshgrid(x,y)
    #fig=pl.figure()
    #ax=Axes3D(fig)
    #ax.plot_surface(X,Y,ob)

    #pl.show()

