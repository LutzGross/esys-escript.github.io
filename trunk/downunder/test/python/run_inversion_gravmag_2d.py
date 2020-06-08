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

"""2D magnetic/gravity joint inversion example using synthetic data"""

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import os
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.weipa import saveSilo

try:
    import esys.ripley
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False
    
try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

@unittest.skipIf(not HAVE_RIPLEY, "Ripley module not available")
class Test_JoinInversion(unittest.TestCase):
    def test_2D_inversion(self):
        logging.getLogger('inv.minimize.MinimizerLBFGS').setLevel(logging.CRITICAL)
        logging.getLogger('inv.JointGravityMagneticInversion').setLevel(logging.CRITICAL)
        # interesting parameters:
        depth_offset = 10 * U.km
        n_humps_h = 5
        n_humps_v = 2
        n_cells_in_data = 30
        full_knowledge = False
        B_b = [2201.*U.Nano*U.Tesla, 31232.*U.Nano*U.Tesla, -41405.*U.Nano*U.Tesla]
        #
        DIM = 2
        n_cells_in_data = max(n_humps_h*7, n_cells_in_data)
        l_data = 100. * U.km
        l_pad = 40. * U.km
        THICKNESS = 20. * U.km
        l_air = 20. * U.km
        n_cells_v = max(
                int((2*l_air+THICKNESS+depth_offset)/l_data*n_cells_in_data + 0.5), 25)


        grav_data=SyntheticData(DataSource.GRAVITY, n_length=n_humps_h, n_depth=n_humps_v,
                depth=THICKNESS+depth_offset, depth_offset=depth_offset,
                DIM=DIM, number_of_elements=n_cells_in_data, length=l_data,
                data_offset=0, full_knowledge=full_knowledge)

        mag_data=SyntheticData(DataSource.MAGNETIC, n_length=n_humps_h, n_depth=n_humps_v,
                depth=THICKNESS+depth_offset, depth_offset=depth_offset,
                DIM=DIM, number_of_elements=n_cells_in_data, length=l_data, B_b=B_b,
                data_offset=0, full_knowledge=full_knowledge, s=l_data/n_humps_h*0.1)

        domainbuilder=DomainBuilder(dim=DIM)
        domainbuilder.addSource(grav_data)
        domainbuilder.addSource(mag_data)
        domainbuilder.setVerticalExtents(depth=l_air+THICKNESS+depth_offset,
                                         air_layer=l_air, num_cells=n_cells_v)
        domainbuilder.setPadding(l_pad)
        domainbuilder.fixDensityBelow(depth=THICKNESS+depth_offset)
        domainbuilder.fixSusceptibilityBelow(depth=THICKNESS+depth_offset)
        domainbuilder.setBackgroundMagneticFluxDensity(B_b)

        inv=JointGravityMagneticInversion()
        inv.setSolverTolerance(1e-3)
        inv.setSolverMaxIterations(500)
        inv.setup(domainbuilder)

        inv.getCostFunction().setTradeOffFactorsModels([1., 0.01])
        inv.getCostFunction().setTradeOffFactorsRegularization(mu = [1.e-2,1.e-2], mu_c=1000.)

        rho_ref = grav_data.getReferenceProperty()
        k_ref = mag_data.getReferenceProperty()

        rho_new, k_new = inv.run()

#        print("rho_new = %s"%rho_new)
#        print("rho = %s"%rho_ref)
#        print("k_new = %s"%k_new)
#        print("k = %s"%k_ref)

        g, chi = inv.getCostFunction().getForwardModel(inv.DENSITY).getSurvey(0)
        B, chi = inv.getCostFunction().getForwardModel(inv.SUSCEPTIBILITY).getSurvey(0)

#        saveSilo(os.path.join(WORKDIR, 'results_joint_2d'),
#                 density=rho_new, density_ref=rho_ref,
#                 susceptibility=k_new, susceptibility_ref=k_ref,
#                 g_data=g, B_data=B, chi=chi)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
