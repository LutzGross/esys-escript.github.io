
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

"""2D gravity inversion example using synthetic data"""

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import os
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.weipa import saveSilo

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

import logging
logger=logging.getLogger('inv')
logger.setLevel(logging.INFO)

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

