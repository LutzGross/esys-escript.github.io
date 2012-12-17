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

from datetime import datetime
import numpy as np
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
MISSING=-9999

# Units of the data array
UNITS='mGal'

# The actual data array, must have shape (NY, NX)
DATA=10*np.random.normal(size=(NY, NX))

##############################################################################
###################### Keep everything below this line #######################
##############################################################################

# conventions used in the file
conventions="CF-1.0, COARDS, Unidata Dataset Discovery v1.0"
# projection string
esri_pe_string="GEOGCS[\\\"GCS_WGS_1984\\\",DATUM[\\\"D_WGS_1984\\\",SPHEROID[\\\"WGS_1984\\\",6378137.0,298.257223563]],PRIMEM[\\\"Greenwich\\\",0.0],UNIT[\\\"Degree\\\",0.0174532925199433]]"
# file log
history=datetime.now().strftime("%d-%m-%Y")+" created using python script"
# license
license="Free to use"

longitude=np.linspace(ORIGIN_X, ORIGIN_X+NX*DELTA_X, NX, endpoint=False)
latitude=np.linspace(ORIGIN_Y, ORIGIN_Y-NY*DELTA_Y, NY, endpoint=False)

o=netcdf_file(FILENAME,'w')
o.cdm_data_type="Grid"
o.Conventions=conventions
o.history=history
o.license=license
o.Metadata_Conventions=conventions
o.summary=SUMMARY
o.title=TITLE
o.createDimension("latitude", NY)
o.createDimension("longitude", NX)

v=o.createVariable("latitude", latitude.dtype, ["latitude"])
v.axis = "Y"
v.long_name="Latitude"
v.standard_name="latitude"
v.units="degrees_north"
v.data[:]=latitude

v=o.createVariable("longitude", longitude.dtype, ["longitude"])
v.axis = "X"
v.long_name="Longitude"
v.standard_name="longitude"
v.units="degrees_east"
v.data[:]=longitude

v=o.createVariable("data", DATA.dtype, ["latitude","longitude"])
v.coordinates="lon lat"
v.esri_pe_string=esri_pe_string
v.long_name="Bouguer_anomaly"
v.missing_value=MISSING
v.units=UNITS
v.data[:]=DATA

o.close()

