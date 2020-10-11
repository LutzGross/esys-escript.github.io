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

"""3D magnetic inversion example using netCDF data"""
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
from esys.escript import saveDataCSV, hasFeature

# Set parameters
DATASET = 'data/MagneticSmall.nc'
DATA_UNITS = U.Nano * U.Tesla
PAD_X = 0.2
PAD_Y = 0.2
thickness = 40. * U.km
l_air = 6. * U.km
n_cells_v = 25
MU = 0.1
# background magnetic field components (B_East, B_North, B_Vertical)
B_b = [2201.*U.Nano*U.Tesla, 31232.*U.Nano*U.Tesla, -41405.*U.Nano*U.Tesla]

COORDINATES=CartesianReferenceSystem()
#COORDINATES=WGS84ReferenceSystem()


def work():
  # Setup and run the inversion
  source=NetCdfData(NetCdfData.MAGNETIC, DATASET, scale_factor=DATA_UNITS, reference_system=COORDINATES)
  db=DomainBuilder(dim=3, reference_system=COORDINATES)
  db.addSource(source)
  db.setVerticalExtents(depth=thickness, air_layer=l_air, num_cells=n_cells_v)
  db.setFractionalPadding(pad_x=PAD_X, pad_y=PAD_Y)
  db.setBackgroundMagneticFluxDensity(B_b)
  db.fixSusceptibilityBelow(depth=thickness)

  inv=MagneticInversion(self_demagnetization=True)
  inv.setSolverTolerance(1e-4)
  inv.setSolverMaxIterations(100)
  inv.fixMagneticPotentialAtBottom(False)
  inv.setup(db)
  inv.getCostFunction().setTradeOffFactorsModels(MU)

  susceptibility = inv.run()
  print("susceptibility = %s"%susceptibility)

  B, w =  db.getMagneticSurveys()[0]
  try:
      saveSilo("result_magnetic.silo", susceptibility=susceptibility, magnetic_anomaly=B, magnetic_weight=w)
      print("Results saved in result_magnetic.silo")
  except:
      print("Failed to save result_magnetic.silo. Possibly no Silo support.")

  saveVTK("result_magnetic.vtu", susceptibility=susceptibility, magnetic_anomaly=B, magnetic_weight=w)
  print("Results saved in result_magnetic.vtu")

  saveDataCSV("result_magnetic.csv", susceptibility=susceptibility, x=susceptibility.getFunctionSpace().getX())
  print("Results saved in result_magnetic.csv")

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

if not hasFeature('netcdf'):
    print("This example requires scipy's netcdf support which does not appear to be installed.")
elif not HAVE_RIPLEY:
    print("Ripley module required but not available")
elif not HAVE_PYPROJ:
    print("This example requires pyproj to be installed.")
else:
    work()
