
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import logging
import os
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.weipa import saveSilo


try: 
   WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
   WORKDIR='.'
logger=logging.getLogger('inv')
logger.setLevel(logging.DEBUG)
handler=logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)

# interesting parameter:
depth_offset=0.*U.km
n_humbs_h= 3
n_humbs_v=1
mu=100
n_cells_in_data=100
full_knowledge=False
spherical=False
# 
n_cells_in_data=max(n_humbs_h*7,n_cells_in_data)
l_data = 100 * U.km
l_pad=40*U.km
THICKNESS=20.*U.km
l_air=20*U.km
n_cells_v=max(int((2*l_air+THICKNESS+depth_offset)/l_data*n_cells_in_data + 0.5), 25)

B_b=simpleGeoMagneticFluxDensity(latitude=-28.5)

                     
grav_data=SyntheticData(DataSource.GRAVITY,n_length=n_humbs_h, n_depth=n_humbs_v, depth=THICKNESS+depth_offset, depth_offset=depth_offset,
                     DIM=2, number_of_elements =n_cells_in_data, length=l_data,
                     data_offset=0,full_knowledge=full_knowledge, spherical=spherical)

mag_data=SyntheticData(DataSource.MAGNETIC,n_length=n_humbs_h, n_depth=n_humbs_v, depth=THICKNESS+depth_offset, depth_offset=depth_offset,
                     DIM=2, number_of_elements =n_cells_in_data, length=l_data, B_b=B_b,
                     data_offset=0,full_knowledge=full_knowledge, spherical=spherical)
                   

domainbuilder=DomainBuilder(dim=2)
domainbuilder.addSource(grav_data)
domainbuilder.addSource(mag_data)
domainbuilder.setVerticalExtents(depth=l_air+THICKNESS+depth_offset, air_layer=l_air, num_cells=n_cells_v)
domainbuilder.setPadding(l_pad)
domainbuilder.fixDensityBelow(depth=THICKNESS+depth_offset)
domainbuilder.fixSusceptibilityBelow(depth=THICKNESS+depth_offset)
domainbuilder.setBackgroundMagneticFluxDensity(B_b)


inv=JointGravityMagneticInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(50)
inv.setup(domainbuilder)


inv.getCostFunction().setTradeOffFactorsModels([10., 1.])
inv.getCostFunction().setTradeOffFactorsRegularization(mu = [1.,1.], mu_c=1.)


rho_new, k_new =inv.run()
rho=grav_data.getReferenceProperty()
k = mag_data.getReferenceProperty()
print "rho_new = ",rho_new
print "rho =", rho
print "k_new = ",k_new
print "k =", k

g, chi = inv.getCostFunction().getForwardModel(inv.DENSITY).getSurvey(0)
B, chi = inv.getCostFunction().getForwardModel(inv.SUSCEPTIBILITY).getSurvey(0)
                                                   
saveSilo(os.path.join(WORKDIR, 'jointinv'), 
     density=rho_new, density_ref=rho,
     susceptability=k_new, susceptability_ref=k,
     g_data=g, B_data=B, chi=chi)

