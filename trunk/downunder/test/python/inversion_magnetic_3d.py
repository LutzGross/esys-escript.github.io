
##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
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

"""3D magnetic inversion example using synthetic data"""

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
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

# interesting parameters:
depth_offset = 10. * U.km
n_humps_h = 5
n_humps_v = 2
mu = 0.1
n_cells_in_data = 50
latitude = -28.5
full_knowledge = False
#
DIM = 3
n_cells_in_data = max(n_humps_h*7, n_cells_in_data)
l_data = 100. * U.km
l_pad = 40. * U.km
THICKNESS = 20. * U.km
l_air = 20. * U.km
n_cells_v = max(
        int((2*l_air+THICKNESS+depth_offset)/l_data*n_cells_in_data + 0.5), 25)

B_b=simpleGeoMagneticFluxDensity(latitude=latitude)

source=SyntheticData(DataSource.MAGNETIC, n_length=n_humps_h, n_depth=n_humps_v,
        depth=THICKNESS+depth_offset, depth_offset=depth_offset,
        DIM=DIM, number_of_elements=n_cells_in_data, length=l_data, B_b=B_b,
        data_offset=0, full_knowledge=full_knowledge, spherical=False)

domainbuilder=DomainBuilder(dim=DIM)
domainbuilder.addSource(source)
domainbuilder.setVerticalExtents(depth=l_air+THICKNESS+depth_offset,
                                 air_layer=l_air, num_cells=n_cells_v)
domainbuilder.setBackgroundMagneticFluxDensity(B_b)
domainbuilder.setPadding(pad_x=l_pad, pad_y=l_pad)
domainbuilder.fixSusceptibilityBelow(depth=THICKNESS+depth_offset)

inv=MagneticInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(50)
inv.setup(domainbuilder)
inv.getCostFunction().setTradeOffFactorsModels(mu)

k_new=inv.run()
k_ref=source.getReferenceProperty()
print("k_new = %s"%k_new)
print("k = %s"%k_ref)
B, chi = inv.getCostFunction().getForwardModel().getSurvey(0)
saveSilo(os.path.join(WORKDIR, 'results_magnetic_3d'), k=k_new, k_ref=k_ref, B=B, chi=chi)

