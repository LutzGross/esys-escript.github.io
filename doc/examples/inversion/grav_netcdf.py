##############################################################################
#
# Copyright (c) 2009-2018 by The University of Queensland
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

"""3D gravity inversion example using netCDF data"""
from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# Import required modules
from esys.downunder import *
from esys.weipa import *
from esys.escript import unitsSI as U
from esys.escript import saveDataCSV

# Set parameters
DATASET = 'data/GravitySmall.nc'
DATA_UNITS = 1e-6 * U.m/(U.sec**2)
PAD_X = 0.2
PAD_Y = 0.2
thickness = 40. * U.km
l_air = 6. * U.km
n_cells_v = 15
MU = 0.1

COORDINATES=CartesianReferenceSystem()
#COORDINATES=WGS84ReferenceSystem()

# Setup and run the inversion
def work():
  source=NetCdfData(NetCdfData.GRAVITY, DATASET, scale_factor=DATA_UNITS, reference_system=COORDINATES)
  db=DomainBuilder(dim=3, reference_system=COORDINATES)
  db.addSource(source)
  db.setVerticalExtents(depth=thickness, air_layer=l_air, num_cells=n_cells_v)
  db.setFractionalPadding(pad_x=PAD_X, pad_y=PAD_Y)
  db.fixDensityBelow(depth=thickness)

  inv=GravityInversion()
  inv.setSolverTolerance(1e-4)
  inv.setSolverMaxIterations(50)
  inv.setup(db)
  inv.getCostFunction().setTradeOffFactorsModels(MU)

  density = inv.run()
  print("density = %s"%density)

  g, w =  db.getGravitySurveys()[0]
  saveVoxet("result.vo", density=density)
  if saveSilo("result_gravity.silo", density=density, gravity_anomaly=g, gravity_weight=w):
      print("Results saved in result_gravity.silo")
  else:
      print("Failed to save result_gravity.silo. Possibly no Silo support.")

  saveVTK("result_gravity.vtu", density=density, gravity_anomaly=g, gravity_weight=w)
  print("Results saved in result_gravity.vtu")

  saveDataCSV("result_gravity.csv", density=density, x=density.getFunctionSpace().getX())
  print("Results saved in result_gravity.csv")
  print("All done. Have a nice day!")

try:
    import pyproj
    HAVE_PYPROJ = True
except ImportError:
    HAVE_PYPROJ = False

try:
    import esys.ripley
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False

if 'NetCdfData' not in dir():
    print("This example requires scipy's netcdf support which does not appear to be installed.")
elif not HAVE_RIPLEY:
    print("Ripley module not available")
elif not HAVE_PYPROJ:
    print("This example requires pyproj to be installed.")
else:
    work()

