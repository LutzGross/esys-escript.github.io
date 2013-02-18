
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

"""3D gravity/magnetic joint inversion example using netCDF data"""

__copyright__="""Copyright (c) 2009-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# Import required modules
from esys.downunder import *
from esys.downunder.minimizers import MinimizerMaxIterReached
from esys.escript import unitsSI as U
from esys.weipa import saveSilo

# Set parameters
MAGNETIC_DATASET = 'magnetic_anomaly.nc'
GRAVITY_DATASET = 'gravity_anomaly.nc'
latitude = -28.5
PAD_X = 0.2
PAD_Y = 0.2
thickness = 40. * U.km
l_air = 6. * U.km
n_cells_v = 25
mu_gravity = 10.
mu_magnetic = 0.1

# Setup and run the inversion
B_b=simpleGeoMagneticFluxDensity(latitude=latitude)
grav_source=NetCdfData(NetCdfData.GRAVITY, GRAVITY_DATASET)
mag_source=NetCdfData(NetCdfData.MAGNETIC, MAGNETIC_DATASET)
db=DomainBuilder(dim=3)
db.addSource(grav_source)
db.addSource(mag_source)
db.setVerticalExtents(depth=thickness, air_layer=l_air, num_cells=n_cells_v)
db.setFractionalPadding(pad_x=PAD_X, pad_y=PAD_Y)
db.setBackgroundMagneticFluxDensity(B_b)
db.fixDensityBelow(depth=thickness)
db.fixSusceptibilityBelow(depth=thickness)

inv=JointGravityMagneticInversion()
inv.setSolverTolerance(1e-4)
inv.setSolverMaxIterations(50)
inv.setup(db)
inv.getCostFunction().setTradeOffFactorsModels([mu_gravity, mu_magnetic])
inv.getCostFunction().setTradeOffFactorsRegularization(mu = [1.,1.], mu_c=1.)

print("Starting inversion, please stand by...")
try:
    density, susceptibility = inv.run()
except MinimizerMaxIterReached as e:
    print(e)
    density, susceptibility = inv.p

print("density = %s"%density)
print("susceptibility = %s"%susceptibility)

# Save results
g, wg =  inv.getCostFunction().getForwardModel(0).getSurvey()
B, wB =  inv.getCostFunction().getForwardModel(1).getSurvey()
saveSilo("result_gravmag.silo", density=density, gravity_anomaly=g[2], gravity_weight=wg[2], susceptibility=susceptibility, magnetic_anomaly=B[2], magnetic_weight=wB[2])
print("Results saved in result_gravmag.silo. Good bye.")

