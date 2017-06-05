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
# example10a.py
# Model of gravitational Potential for an infinite strike horizontal
# cylinder. This example tests the boundary condition criterion.

#######################################################EXTERNAL MODULES
# To solve the problem it is necessary to import the modules we require.
import matplotlib
matplotlib.use('agg') #It's just here for automated testing

from esys.escript import * # This imports everything from the escript library
from esys.escript.unitsSI import * 
from esys.escript.linearPDEs import LinearPDE # This defines LinearPDE as LinearPDE
from esys.weipa import saveVTK # This imports the VTK file saver from weipa
import os, sys #This package is necessary to handle saving our data.
from math import pi, sqrt, sin, cos

from esys.escript.pdetools import Projector, Locator
from cblib import toRegGrid
from esys.finley import ReadGmsh

import pylab as pl #Plotting package
import numpy as np

try:
    # This imports the rectangle domain function 
    from esys.finley import Rectangle
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
    #PDE related
    rho=200.0                   #the density contrast of the source
    G=6.67300*10E-11            #gravitational constant
    R=100.                      #radius of the source body
    ################################################ESTABLISHING PARAMETERS
    #the folder to put our outputs in, leave blank "" for script path 
    save_path= os.path.join("data","example10")
    mesh_path = "data/example10m"
    #ensure the dir exists
    mkDir(save_path)

    ####################################################DOMAIN CONSTRUCTION
    #Geometric and material property related variables.
    domain=ReadGmsh(os.path.join(mesh_path,'example10m.msh'),2) # create the domain
    x=domain.getX()
    #Domain related dimensions from the imported file.
    mx = Lsup(x)*m #meters - model length
    my = Lsup(x)*m #meters - model width
    rholoc=[mx/2,my/2]

    mask=wherePositive(R-length(x-rholoc))
    rhoe=rho*mask
    kro=kronecker(domain)

    #Define the analytic solution.
    def analytic_gz(x,z,R,drho):
        G=6.67300*10E-11
        return G*2*np.pi*R*R*drho*(z/(x*x+z*z))

    #Define the boundary conditions.
    q=whereZero(x[1]-my)+whereZero(x[1])+whereZero(x[0])+whereZero(x[0]-mx)
    ###############################################ESCRIPT PDE CONSTRUCTION
    mypde=LinearPDE(domain)
    mypde.setValue(A=kro,Y=4.*3.1415*G*rhoe,q=q,r=0.)
    mypde.setSymmetryOn()
    sol=mypde.getSolution()

    g_field=grad(sol) #The gravitational acceleration g.
    g_fieldz=g_field*[0,1] #The vertical component of the g field.
    gz=length(g_fieldz) #The magnitude of the vertical component.
    # Save the output to file.
    saveVTK(os.path.join(save_path,"ex10d.vtu"),\
            grav_pot=sol,g_field=g_field,g_fieldz=g_fieldz,gz=gz,rho=rhoe)

    ################################################MODEL SIZE SAMPLING
    smoother = Projector(domain)  #Function smoother.
    z=1000.                       #Distance of profile from source.
    sol_angz=[]                   #Array for analytic gz
    sol_anx=[]                    #Array for x
    #calculate analytic gz and x location.
    for x in range(int(-mx/2.),int(mx/2.),10):
        sol_angz.append(analytic_gz(x,z,R,rho))
        sol_anx.append(x+mx/2)
    #save analytic solution
    pl.savetxt(os.path.join(save_path,"ex10d_as.asc"),sol_angz)

    sol_escgz=[]                  #Array for escript solution for gz
    #calculate the location of the profile in the domain
    for i in range(0,len(sol_anx)):
        sol_escgz.append([sol_anx[i],my/2.-z])

    rec=Locator(gz.getFunctionSpace(),sol_escgz) #location to record
    psol=rec.getValue(smoother(gz))              #extract the solution
    #save X and escript solution profile to file
    pl.savetxt(os.path.join(save_path,"ex10d_%05d.asc"%mx),psol)
    pl.savetxt(os.path.join(save_path,"ex10d_x%05d.asc"%mx),sol_anx)

    #plot the pofile and analytic solution using example10q.py

