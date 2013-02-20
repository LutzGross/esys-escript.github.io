##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
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

"""This example show how to display netCDF input data with matplotlib"""

from matplotlib import pyplot as plt
import numpy as np
import sys
from scipy.io import netcdf_file

# input filename
if len(sys.argv)>1:
    FILENAME=sys.argv[1]
else:
    FILENAME='data/QLDWest_grav.nc'

f=netcdf_file(FILENAME, 'r')
NY=f.dimensions["latitude"]
NX=f.dimensions["longitude"]

latitude=f.variables["latitude"]
y_label=latitude.long_name
y_units=latitude.units
latitude=latitude[:]

longitude=f.variables["longitude"]
x_label=longitude.long_name
x_units=longitude.units
longitude=longitude[:]

DATA=f.variables["onshore_only_Bouguer_geodetic"]
data_label=DATA.long_name
UNITS=DATA.units
DATA=DATA[:]

f.close()

# add one more point in each dimension (because we fill cells):
ll=2*longitude[-1]-longitude[-2]
longitude=np.resize(longitude, len(longitude)+1)
longitude[-1]=ll
ll=2*latitude[-1]-latitude[-2]
latitude=np.resize(latitude, len(latitude)+1)
latitude[-1]=ll

lx=longitude[-1]-longitude[0]
ly=latitude[-1]-latitude[0]
x,y=np.meshgrid(longitude, latitude)
plt.figure(figsize=(6*lx/ly+1, 6), dpi=100)
plt.pcolor(x, y, DATA)
locs,_=plt.xticks()
plt.xticks(locs, map(lambda x:"%g"%x, locs))
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.axis('tight')
plt.title("%s (%s)"%(data_label,UNITS))
plt.colorbar()
plt.show()

