
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

"""this a very simple test for inversion of acoustic data in the freqency domain"""

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
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
N_OMEGA=3
# these are the taeget P-velocity and Q index
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
SIGMA=1./( V_P * complex(1., -1./(2.*Q)))**2
S={}
for o in omegas:
    S[o]=complex(sqrt(o), sin(o))
    
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
reg_mask=Data(0.,(2,), Solution(dom))
reg_mask[0] = whereZero(z+inf(z))        # mask for locations where velocity and Q index are assumed to be known
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

m=AcousticVelocityMapping(V_P_prior, Q_prior)
# finally we can set up the cost function:
cf=InversionCostFunction(regularization,
             (m,) , tuple(forwardmodels) )
             
# and then start running the show:
solver=MinimizerLBFGS()
solver.setCostFunction(cf)
solver.setTolerance(1e-4)
solver.setMaxIterations(50)
solver.run(Data(0.,(2,),Solution(dom)))
m=solver.getResult()
sigma = cf.getProperties(m)
  
print sigma




1/0
#=============================================================
# interesting parameters:
n_humps_h = 3
n_humps_v = 1
mu = 100
n_cells_in_data = 100
# ignore:
full_knowledge = False
depth_offset = 0. * U.km
#
DIM = 2
n_cells_in_data = max(n_humps_h*7, n_cells_in_data)
l_data = 100. * U.km
l_pad = 40. * U.km
THICKNESS = 20. * U.km
l_air = 20. * U.km
n_cells_v = max(
        int((2*l_air+THICKNESS+depth_offset)/l_data*n_cells_in_data + 0.5), 25)


source=SyntheticData(DataSource.GRAVITY, n_length=n_humps_h, n_depth=n_humps_v,
        depth=THICKNESS+depth_offset, depth_offset=depth_offset,
        DIM=DIM, number_of_elements=n_cells_in_data, length=l_data,
        data_offset=0, full_knowledge=full_knowledge)

domainbuilder=DomainBuilder(dim=DIM)
domainbuilder.addSource(source)
domainbuilder.setVerticalExtents(depth=l_air+THICKNESS+depth_offset,
                                 air_layer=l_air, num_cells=n_cells_v)
domainbuilder.setPadding(l_pad)
domainbuilder.fixDensityBelow(depth=THICKNESS+depth_offset)

inv=GravityInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(50)
inv.setup(domainbuilder)
inv.getCostFunction().setTradeOffFactorsModels(mu)


rho_new = inv.run()
rho_ref = source.getReferenceProperty()
print("rho_new = %s"%rho_new)
print("rho = %s"%rho_ref)
g, chi = inv.getCostFunction().getForwardModel().getSurvey(0)
saveSilo(os.path.join(WORKDIR, 'results_gravity_2d'), density=rho_new, density_ref=rho_ref, g=g, chi=chi)

