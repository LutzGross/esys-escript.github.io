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

"""This example shows how to display netCDF input data with matplotlib"""

from __future__ import division, print_function
import matplotlib
# The following line is here to allow automated testing. Remove or comment if
# you would like to display the final plot in a window instead.
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import sys
try:
    from scipy.io import netcdf_file

    # input filename
    if len(sys.argv)>1:
        FILENAME=sys.argv[1]
    else:
        FILENAME='data/QLDWestMagnetic.nc'

    f=netcdf_file(FILENAME, 'r')
    for latvar in "latitude", "lat", "y":
        try:
            NY=f.dimensions[latvar]
            break
        except KeyError:
            pass
    for lonvar in "longitude", "lon", "x":
        try:
            NX=f.dimensions[lonvar]
            break
        except KeyError:
            pass

    latitude=f.variables[latvar]
    y_label=latitude.long_name
    y_units=latitude.units
    latitude=latitude[:]

    longitude=f.variables[lonvar]
    x_label=longitude.long_name
    x_units=longitude.units
    longitude=longitude[:]

    DATA=f.variables["magmap_V5_2010"]
    data_label=DATA.long_name
    try:
        UNITS=DATA.units
    except AttributeError:
        UNITS="nT"

    # convert NaN (not-a-number) values to -1000 so the plotting works
    DATA=np.where(np.isnan(DATA[:]), -1000, DATA[:])

    # make sure data is in native format
    if DATA.dtype.byteorder != '=':
        DATA=np.array(DATA, dtype=np.float32)

    # add one more point in each dimension (because we fill cells):
    ll=2*longitude[-1]-longitude[-2]
    longitude=np.resize(longitude, len(longitude)+1)
    longitude[-1]=ll
    ll=2*latitude[-1]-latitude[-2]
    latitude=np.resize(latitude, len(latitude)+1)
    latitude[-1]=ll

    lx=abs(longitude[-1]-longitude[0])
    ly=abs(latitude[-1]-latitude[0])
    x,y=np.meshgrid(longitude, latitude)
    plt.figure(figsize=(6*lx/ly+1, 6), dpi=100)
    plt.pcolor(x, y, DATA)
    locs,_=plt.xticks()
    plt.xticks(locs, list(map(lambda x:"%g"%x, locs)))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.axis('tight')
    plt.title("%s (%s)"%(data_label,UNITS))
    plt.colorbar()

    plt.show()
    plt.savefig("netcdf_plot.png")
    # Now we can close the file
    f.close()

except ImportError:
    print("The scipy module was not found but is required to read netCDF files. Exiting...")

