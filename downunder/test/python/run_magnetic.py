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

"""2D magnetic inversion example using synthetic data"""

from __future__ import print_function, division

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
class Test_MagneticInversion2D(unittest.TestCase):
    def test_2D_inversion(self):
        logging.getLogger('inv.minimizer.MinimizerLBFGS').setLevel(logging.CRITICAL)
        #logging.getLogger('inv.MagneticInversion').setLevel(logging.CRITICAL)
        logging.getLogger('inv').setLevel(logging.CRITICAL)
        # interesting parameters:
        depth_offset = 0. * U.km
        n_humps_h = 1
        n_humps_v = 1
        mu = 1.
        n_cells_in_data = 200
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


        source=SyntheticData(DataSource.MAGNETIC, n_length=n_humps_h, n_depth=n_humps_v,
                depth=THICKNESS+depth_offset, depth_offset=depth_offset,
                DIM=DIM, number_of_elements=n_cells_in_data, length=l_data, B_b=B_b,
                data_offset=0, full_knowledge=full_knowledge)

        domainbuilder=DomainBuilder(dim=DIM)
        domainbuilder.addSource(source)
        domainbuilder.setVerticalExtents(depth=l_air+THICKNESS+depth_offset,
                                         air_layer=l_air, num_cells=n_cells_v)
        domainbuilder.setBackgroundMagneticFluxDensity(B_b)
        domainbuilder.setPadding(l_pad)
        domainbuilder.fixSusceptibilityBelow(depth=THICKNESS+depth_offset)

        inv=MagneticInversion()
        inv.setSolverTolerance(1e-4)
        inv.setSolverMaxIterations(70)
        inv.setup(domainbuilder)
        inv.getCostFunction().setTradeOffFactorsModels(mu)

        k_new=inv.run()
        k_ref=source.getReferenceProperty()
        B, chi = inv.getCostFunction().getForwardModel().getSurvey(0)
        #saveSilo(os.path.join(WORKDIR, 'results_magnetic_2d'), k=k_new, k_ref=k_ref, B=B, chi=chi)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
