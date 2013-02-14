
##############################################################################
#
# Copyright (c) 2009-2013 by University of Queensland
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

"""3D gravity inversion example using netCDF data"""

__copyright__="""Copyright (c) 2009-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# Import required modules
from esys.downunder import *
from esys.weipa import *
from esys.escript import unitsSI as U

# Set parameters
DATASET = 'data/QLDWest_grav.nc'
PAD_X = 0.2
PAD_Y = 0.2
thickness = 40. * U.km
l_air = 6. * U.km
n_cells_v = 25
MU = 0.1

# Setup and run the inversion
source=NetCdfData(NetCdfData.GRAVITY, DATASET)
db=DomainBuilder(dim=3)
db.addSource(source)
db.setVerticalExtents(depth=thickness, air_layer=l_air, num_cells=n_cells_v)
db.setFractionalPadding(pad_x=PAD_X, pad_y=PAD_Y)
db.fixDensityBelow(depth=thickness)

inv=GravityInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(50)
inv.setup(db)
inv.getCostFunction().setTradeOffFactorsModels(MU)

print("Starting inversion, please stand by...")
rho = inv.run()
print("density = %s"%density)

print("Results saved in result0.silo:")
g, w =  inv.getCostFunction().getForwardModel().getSurvey(0)
saveSilo("result0.silo", density=rho, gravity_anomaly=g[2], gravity_weight=w[2])

print("Results saved in result0.vtu:")
saveSilo("result0.vtu", density=rho)

print("Results saved in result0.csv:")
saveCSV("result0.csv", density=rho, x=rho.getFunctionSpace().getX())

print("All done. Have a nice day.!")
