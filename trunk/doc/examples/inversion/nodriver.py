##############################################################################
#
# Copyright (c) 2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.downunder import *
import esys.escript.unitsSI as U
from esys.weipa import saveSilo

import numpy as np


# Set parameters

DATASET_GRAV='bouguer_anomaly.nc'
DATASET_MAG='bouguer_anomaly.nc'
latitude=-28.5

thickness = 40. * U.km # below surface
l_air = 6. * U.km      # above surface
n_cells_v = 25         # number of cells in vertical direction   

# apply 20% padding
PAD_X = 0.2
PAD_Y = 0.2

mu = 10.

# read data:
source_g=NetCdfData(NetCdfData.GRAVITY, DATASET_GRAV)
source_m=NetCdfData(NetCdfData.MAGNETIC, DATASET_MAG)
B_b=simpleGeoMagneticFluxDensity(latitude=latitude)


# create domain:
db=DomainBuilder(dim=3)
db.addSource(source_g)
db.addSource(source_m)
db.setVerticalExtents(depth=thickness, air_layer=l_air, num_cells=n_cells_v)
db.setFractionalPadding(pad_x=PAD_X, pad_y=PAD_Y)
db.fixDensityBelow(depth=thickness)

dom=db.getDomain()
DIM=dom.getDim()

# create mappings with standard parameters
rho_mapping=DensityMapping(dom)
k_mapping=SusceptibilityMapping(dom)
         
# create regularization with two level set functions: 
reg_mask=Data(0.,(2,), Solution(dom))   
reg_mask[0] = domainbuilder.getSetDensityMask()          # mask for locations where m[0]~rho is known
reg_mask[1] = domainbuilder.getSetSusceptibilityMask()   # mask for locations where m[0]~k is known
regularization=Regularization(dom, numLevelSets=2, 
                               w1=np.ones((2,DIM)),      # consider gradient terms
                               wc=[[0,0],[1,0]],         # and cross gradient term
                               location_of_set_m=reg_mask)
                               
# create forward model for gravity
# get data with deviation
g,sigma_g=domainbuilder.getGravitySurveys()[0]
# turn the scalars into vectors (vertical direction)
d=kronecker(DIM)[DIM-1]
w=safeDiv(1., sigma_g)

gravity_model=GravityModel(dom, w*d, g*d)
gravity_model.rescaleWeights(rho_scale=rho_mapping.getTypicalDerivative())

# create forward model for magnetic
d=normalize(B_b) # direction of measurement 
B,sigma_B=domainbuilder.getMagneticSurveys()[0]
w=safeDiv(1., sigma_B)

magnetic_model=MagneticModel(dom, w*d, B*d, B_b)
magnetic_model.rescaleWeights(k_scale=k_mapping.getTypicalDerivative())


# finally we can set up the cost function:
cf=InversionCostFunction(regularization, 
             ((rho_mapping,0), (k_mapping, 1)),
             ((gravity_model,0), (magnetic_model,1)) )

cf.setTradeOffFactorsModels(mu)

# sun solver:
solver=MinimizerLBFGS()
solver.setSolverTolerance(1e-4)
solver.setSolverMaxIterations(50)
m=solver.run(Data(0.,(2,),Solution(dom)))
density, susceptibility =getProperties(m)


# write everything to file:
saveSilo("results.silo",
         density=density, susceptability=susceptibility,
         g_data=g, sigma_g=g_sigma_g, B_data=B,  sigma_B=sigma_B)

