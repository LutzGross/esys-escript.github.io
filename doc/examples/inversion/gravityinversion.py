
##############################################################################
#
# Copyright (c) 2009-2012 by University of Queensland
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

__copyright__="""Copyright (c) 2009-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# Import required modules
from esys.downunder import NetCdfData, DomainBuilder, GravityInversion
from esys.escript import unitsSI as U
from esys.weipa import saveSilo

# Set parameters
PAD_X=0.2
PAD_Y=0.2
DATASET='bouguer_anomaly.nc'
DEPTH=40 * U.km
AIR=6 * U.km
NE_Z=25

# Setup and run the inversion
source=NetCdfData(NetCdfData.GRAVITY, DATASET)
db=DomainBuilder()
db.addSource(source)
db.setVerticalExtents(depth=DEPTH, air_layer=AIR, num_cells=NE_Z)
db.setFractionalPadding(PAD_X, PAD_Y)
inv=GravityInversion()
inv.setup(db)
g, w = inv.getForwardModel().getSurvey(0)
density=inv.run()

# Save results
saveSilo("result.silo", density=density, gravity_anomaly=g[2], gravity_weight=w[2])

