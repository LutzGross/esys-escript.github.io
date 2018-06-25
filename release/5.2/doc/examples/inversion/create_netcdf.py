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

"""
This example shows how to create a netCDF input file that is suitable for
inversions in esys.downunder. 
"""
from __future__ import division, print_function

import sys
from datetime import datetime
import numpy as np
try:
    from scipy.io import netcdf_file

    # output filename
    FILENAME='output.nc'

    # a concise title and summary of the dataset
    TITLE="custom_data"
    SUMMARY="Bouguer gravity anomaly data"

    # Origin longitude (degrees east) and latitude (degrees north)
    ORIGIN_X=130.2
    ORIGIN_Y=-29.1

    # spacing in longitude,latitude direction (degrees)
    DELTA_X=0.05
    DELTA_Y=0.05

    # Number of data points in longitude,latitude direction
    NX=20
    NY=10

    # Dummy value (for unset areas)
    MISSING=np.nan

    # Data error (can be constant or variable over the data points)
    SIGMA = 3.

    # The actual data array, must have shape (NY, NX).
    # These are just some random numbers.
    DATA = 10*np.random.normal(size=(NY, NX), scale=SIGMA)
    ERROR = np.ones((NY, NX)) * SIGMA

    # Use MISSING where data is invalid/not available, example:
    DATA[NY//2,NX//2]=MISSING
    ERROR[NY//2,NX//2]=MISSING

    ##############################################################################
    ###################### Keep everything below this line #######################
    ##############################################################################
    # conventions used in the file
    conventions="CF-1.0, COARDS, Unidata Dataset Discovery v1.0"
    # file log
    history=datetime.now().strftime("%d-%m-%Y")+" created using python script"
    # license
    license="Free to use"

    # Create the output file and write a few metadata entries
    o=netcdf_file(FILENAME,'w')
    o.Conventions=conventions
    o.Metadata_Conventions=conventions
    o.history=history
    o.license=license
    o.summary=SUMMARY
    o.title=TITLE

    # Create longitude dimension and variable
    longitude=np.linspace(ORIGIN_X, ORIGIN_X+NX*DELTA_X, NX, endpoint=False)
    o.createDimension("longitude", NX)
    v=o.createVariable("longitude", longitude.dtype, ["longitude"])
    v.data[:]=longitude
    v.units="degrees_east"
    v.long_name="Longitude"

    # Create latitude dimension and variable
    latitude=np.linspace(ORIGIN_Y-NY*DELTA_Y, ORIGIN_Y, NY, endpoint=False)
    o.createDimension("latitude", NY)
    v=o.createVariable("latitude", latitude.dtype, ["latitude"])
    v.data[:]=latitude
    v.units="degrees_north"
    v.long_name="Latitude"

    # Create the main data variable
    v=o.createVariable("data", DATA.dtype, ["latitude","longitude"])
    v.missing_value=MISSING
    v.data[:]=DATA
    v.long_name="Bouguer_anomaly"

    # Create the error variable (can be omitted)
    v=o.createVariable("error", ERROR.dtype, ["latitude","longitude"])
    v.missing_value=MISSING
    v.data[:]=ERROR

    # Close the file
    o.close()

except ImportError:
    print("The scipy module was not found but is required to write netCDF files. Exiting...")

