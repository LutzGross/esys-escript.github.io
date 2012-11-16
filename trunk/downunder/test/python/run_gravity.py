
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

features=[SmoothAnomaly(lx=30*U.km, ly=20*U.km, lz=18.*U.km, \
     x=8*U.km, y=3*U.km, depth=2.5*U.km, v_inner=200., v_outer=1e-6),\
          SmoothAnomaly(lx=25*U.km, ly=20*U.km, lz=20*U.km,
     x=30*U.km, y=1*U.km, depth=18*U.km, v_inner=-200., v_outer=1e-6),\
          SmoothAnomaly(lx=30*U.km, ly=20*U.km, lz=18.*U.km, \
     x=68*U.km, y=3*U.km, depth=5*U.km, v_inner=200., v_outer=1e-6)]

logger=logging.getLogger('inv')
logger.setLevel(logging.DEBUG)
handler=logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
source=SyntheticData(DataSource.GRAVITY, DIM=2, NE=60, l=100*U.km, features=features)
domainbuilder=DomainBuilder(dim=2)
domainbuilder.addSource(source)
domainbuilder.setPadding(10)
domainbuilder.setVerticalExtents(depth=30*U.km, air_layer=10*U.km, num_cells=16)

inv=GravityInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(40)
inv.setSolverOptions(initialHessian=10)
inv.setWeights(mu_model=1e-3)
inv.setup(domainbuilder)

rho_new=inv.run()
g, chi = inv.getForwardModel().getSurvey(0)
saveSilo(os.path.join(WORKDIR, 'gravinv'), density=rho_new, density_ref=source.getReferenceDensity(), g=g, chi=chi)

