
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

from __future__ import print_function, division

"""this a very simple test for inversion of acoustic data in the freqency domain



"""

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import os
from esys.escript import *
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.weipa import saveSilo
import numpy as np

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

import logging
logger=logging.getLogger('inv')
logger.setLevel(logging.INFO)

# frequence range to be used for the inversion
OMEGA_MIN=1*U.Hz
OMEGA_MAX=10*U.Hz
N_OMEGA=1
# these are the target P-velocity and Q index
V_P=2.*U.km
Q=10
# these are the reference values to be used
V_P_prior=1.*U.km
Q_prior=100
#
#   cross gradient weight
#
NU_C=0
#
#   geometry:
#
WIDTH=4.*U.km
NE_X=40.
DEPTH=3*U.km
NE_Z=int(NE_X/WIDTH*DEPTH+0.5)

domainbuilder=DomainBuilder(dim=2)
domainbuilder.setVerticalExtents(depth=DEPTH, air_layer=0., num_cells=NE_Z)

# generate list of omegas:
omegas=np.linspace(OMEGA_MIN, OMEGA_MAX, num=N_OMEGA, endpoint=True)
# create data:
#
#   - laplace u - o**2*SIGMA * u = 1. 
# 
#    -> u=1./o**2/SIGMA = data
#
#    to test the scaling an extra factor S is introduced.
#
SIGMA=1./( V_P * complex(1., -1./(2.*Q)))**2
SIGMA_prior=1./( V_P_prior * complex(1., -1./(2.*Q_prior)))**2
S={}
for o in omegas:
    S[o]=complex(sqrt(o)**2, sin(o))*0+1
    
    DD=-1./(SIGMA*o**2)*S[o]
    data=np.ones((NE_X+1,))*DD
    domainbuilder.addSource(NumpyData(NumpyData.ACOUSTIC, data, error=1., length=WIDTH, tags=[SeismicSource(x=WIDTH/2,omega=o)]))
COORDINATES=None

#====================== THIS IS WERE THE INVERSION STARTS:

tol=1.e-8
saveMemory=True
scaleAllF=False
fixBottomPressure=False
source_width=0

dom=domainbuilder.getDomain()
DIM=dom.getDim()
  
  # create regularization with two level set functions:
z=Solution(dom).getX()[DIM-1]
x=Solution(dom).getX()[0]
reg_mask=Data(0.,(2,), Solution(dom))
reg_mask[0] = whereZero(z-inf(z))        # mask for locations where velocity and Q index are assumed to be known
reg_mask[1] = reg_mask[0]  

regularization=Regularization(dom, numLevelSets=2,
                               w1=np.ones((2,DIM)), # consider gradient terms
                               wc=[[0,NU_C],[0,0]],    # and cross-gradient term
                               coordinates=domainbuilder.getReferenceSystem(),
                               location_of_set_m=reg_mask)

X=Function(dom).getX()
forwardmodels=[]  # this will held the list of forward models 
for ds in domainbuilder.getTags():
    if ds.getPower():
        scaleF=False
        F0=[ds.getPower().real, ds.getPower().imag]
    else:
        scaleF=True
        F0=[1.,0.]

    data, error = domainbuilder.getSurveys(DataSource.ACOUSTIC, [ ds ])[0]
    w=safeDiv(1., error)
    if source_width > 0:
        F=F0*exp(- ((X[DIM-1]-ds.getElevation())/source_width)**2)
        for i in range(DIM-1):
            F*=exp(- ((X[i]-ds.getLocation()[i])/source_width)**2)
    else:
        F=Data(F0, X.getFunctionSpace())
        
    acw=AcousticWaveForm(dom, \
            ds.getFrequency(), w, data, F, \
            coordinates=domainbuilder.getReferenceSystem(), \
            fixAtBottom=fixBottomPressure, \
            tol=tol, \
            saveMemory=saveMemory, \
            scaleF=scaleF)
    logger.debug("power for frequency %s Hz for source at %s added."%(ds.getFrequency, ds.getLocation()))
    forwardmodels.append((acw, 0))

vm=AcousticVelocityMapping(V_P, Q)
# 

z=Solution(dom).getX()[DIM-1]
x=Solution(dom).getX()[0]
z_bar=(z-inf(z))/(sup(z)-inf(z))
x_bar=(x-inf(x))/(sup(x)-inf(x))
s0=((SIGMA_prior.real-SIGMA.real)*z_bar+SIGMA.real)+SIGMA.real*0.1*sin((z_bar*x_bar)*15)
s1=((SIGMA_prior.imag-SIGMA.imag)*z_bar+SIGMA.imag)+SIGMA.imag*0.1*sin((z_bar*x_bar)*15)
m0=vm.getInverse(s0*[1,0]+s1*[0,1])
print "m0=",m0[0]
print "m1=",m0[1]
print "s0=",s0
print "s1=",s1
print "SIGMA=",SIGMA
# finally we can set up the cost function:
cf=InversionCostFunction(regularization,
             (vm,) , tuple(forwardmodels) )
cf.setTradeOffFactorsRegularization(mu=[1e-1,1e-1], mu_c=None)
             
# and then start running the show:
solver=MinimizerLBFGS()
solver.setCostFunction(cf)
solver.setTolerance(1e-6)
solver.setMaxIterations(180)
solver.run(m0)
m=solver.getResult()
sigma = cf.getProperties(m)
  
print "m0=",m[0]
print "m1=",m[1]
print "s0=",sigma[0]
print "s1=",sigma[1]
print "SIGMA=",SIGMA


