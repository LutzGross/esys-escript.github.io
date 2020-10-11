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
from __future__ import division, print_function

"""3D gravity/magnetic joint inversion example using netCDF data"""

__copyright__="""Copyright (c) 2009-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# Import required modules
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.escript import saveDataCSV
from esys.weipa import *

# Set parameters
MAGNETIC_DATASET = 'data/MagneticSmall.nc'
MAG_UNITS = U.Nano * U.Tesla
GRAVITY_DATASET = 'data/GravitySmall.nc'
GRAV_UNITS = 1e-6 * U.m/(U.sec**2)
# background magnetic field components (B_East, B_North, B_Vertical)
B_b = [2201.*U.Nano*U.Tesla, 31232.*U.Nano*U.Tesla, -41405.*U.Nano*U.Tesla]
PAD_X = 0.2
PAD_Y = 0.2
thickness = 40. * U.km
l_air = 6. * U.km
n_cells_v = 25
mu_gravity = 1.
mu_magnetic = 0.1
COORDINATES=CartesianReferenceSystem()
#COORDINATES=WGS84ReferenceSystem()


def work():
  # Setup and run the inversion
  grav_source=NetCdfData(NetCdfData.GRAVITY, GRAVITY_DATASET, scale_factor=GRAV_UNITS, reference_system=COORDINATES)
  mag_source=NetCdfData(NetCdfData.MAGNETIC, MAGNETIC_DATASET, scale_factor=MAG_UNITS, reference_system=COORDINATES)
  db=DomainBuilder(dim=3, reference_system=COORDINATES)
  db.addSource(grav_source)
  db.addSource(mag_source)
  db.setVerticalExtents(depth=thickness, air_layer=l_air, num_cells=n_cells_v)
  db.setFractionalPadding(pad_x=PAD_X, pad_y=PAD_Y)
  db.setBackgroundMagneticFluxDensity(B_b)
  db.fixDensityBelow(depth=thickness)
  db.fixSusceptibilityBelow(depth=thickness)

  inv=StrongJointGravityMagneticInversion()
  inv.setSolverTolerance(1e-4)
  inv.setSolverMaxIterations(50)
  inv.setup(db)
  inv.getCostFunction().setTradeOffFactorsModels([mu_gravity, mu_magnetic])
  inv.getCostFunction().setTradeOffFactorsRegularization(mu = 1.)

  density, susceptibility = inv.run()
  print("density = %s"%density)
  print("susceptibility = %s"%susceptibility)

  g, wg = db.getGravitySurveys()[0]
  B, wB = db.getMagneticSurveys()[0]
  try:
      saveSilo("result_gravmag_strong.silo", density=density, gravity_anomaly=g, gravity_weight=wg, susceptibility=susceptibility, magnetic_anomaly=B,   magnetic_weight=wB)
      print("Results saved in result_gravmag_strong.silo")
  except:
      print("Failed to save silo file. Possibly no Silo support.")

  saveVTK("result_gravmag_strong.vtu", density=density, gravity_anomaly=g, gravity_weight=wg, susceptibility=susceptibility, magnetic_anomaly=B,   magnetic_weight=wB)
  print("Results saved in result_gravmag_strong.vtu")

  print("All done. Have a nice day!")

if 'NetCdfData' in dir():
  work()
else:
  print("This example requires scipy's netcdf support which does not appear to be installed.")

