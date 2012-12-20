
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
     x=8*U.km, y=3*U.km, depth=2.5*U.km, v_inner=2., v_outer=1e-6),\
          SmoothAnomaly(lx=25*U.km, ly=20*U.km, lz=20*U.km,
     x=30*U.km, y=1*U.km, depth=18*U.km, v_inner=-1., v_outer=1e-6),\
          SmoothAnomaly(lx=30*U.km, ly=20*U.km, lz=18.*U.km, \
     x=68*U.km, y=3*U.km, depth=5*U.km, v_inner=20., v_outer=1e-6)]

B_b=simpleGeoMagneticFluxDensity(latitude=-28.5)

logger=logging.getLogger('inv')
logger.setLevel(logging.INFO)
handler=logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)
source=SyntheticFeatureData(DataSource.MAGNETIC, DIM=2, number_of_elements=30, length=100*U.km, features=features, B_b=B_b)

domainbuilder=DomainBuilder(dim=2)
domainbuilder.addSource(source)
domainbuilder.setElementPadding(10)
domainbuilder.setVerticalExtents(depth=30*U.km, air_layer=10*U.km, num_cells=16)
domainbuilder.setBackgroundMagneticFluxDensity(B_b)

inv=MagneticInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(10)
inv.setup(domainbuilder)
k_new = inv.run()
B, chi = inv.getCostFunction().getForwardModel().getSurvey(0)
saveSilo(os.path.join(WORKDIR, 'maginv'), sus=k_new, sus_ref=source.getReferenceProperty(), B=B, chi=chi)

