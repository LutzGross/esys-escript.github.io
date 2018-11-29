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

"""This example shows how to display ER Mapper raster data with matplotlib"""
from __future__ import division, print_function

import matplotlib
# The following line is here to allow automated testing. Remove or comment if
# you would like to display the final plot in a window instead.
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import sys

# input filename
if len(sys.argv)>1:
    FILENAME=sys.argv[1]
else:
    FILENAME='data/QLDWestGravity.ers'


if FILENAME[-4:]=='.ers': FILENAME=FILENAME[:-4]
ersfn=FILENAME+'.ers'
metadata=open(ersfn,'r').readlines()
# parse metadata
start=-1
for i in range(len(metadata)):
    if metadata[i].strip() == 'DatasetHeader Begin':
        start=i+1
if start==-1:
    raise RuntimeError('Invalid ERS file (DatasetHeader not found)')

md_dict={}
section=[]
for i in range(start, len(metadata)):
    line=metadata[i].strip()
    if line[-6:].strip() == 'Begin':
        section.append(line[:-6].strip())
    elif line[-4:].strip() == 'End':
        if len(section)>0:
            section.pop()
    else:
        vals=line.split('=')
        if len(vals)==2:
            key = vals[0].strip()
            value = vals[1].strip()
            fullkey='.'.join(section+[key])
            md_dict[fullkey]=value

try:
    if md_dict['RasterInfo.CellType'] == 'IEEE4ByteReal':
        datatype=np.float32
    elif md_dict['RasterInfo.CellType'] == 'IEEE8ByteReal':
        datatype=np.float64
    elif md_dict['RasterInfo.CellType'] == 'Signed32BitInteger':
        datatype=np.int32
    else:
        raise RuntimeError('Unsupported data type '+md_dict['RasterInfo.CellType'])
except KeyError:
    print("Cell type not specified. Assuming IEEE4ByteReal.")
    datatype=np.float32

try:
    NX = int(md_dict['RasterInfo.NrOfCellsPerLine'])
    NY = int(md_dict['RasterInfo.NrOfLines'])
except:
    raise RuntimeError("Could not determine extents of data")

try:
    spacingX = float(md_dict['RasterInfo.CellInfo.Xdimension'])
    spacingY = float(md_dict['RasterInfo.CellInfo.Ydimension'])
except:
    raise RuntimeError("Could not determine cell dimensions")

try:
    if md_dict['CoordinateSpace.CoordinateType']=='EN':
        originX = float(md_dict['RasterInfo.RegistrationCoord.Eastings'])
        originY = float(md_dict['RasterInfo.RegistrationCoord.Northings'])
        labelX = "Easting"
        labelY = "Northing"
    elif md_dict['CoordinateSpace.CoordinateType']=='LL':
        originX = float(md_dict['RasterInfo.RegistrationCoord.Longitude'])
        originY = float(md_dict['RasterInfo.RegistrationCoord.Latitude'])
        labelX = "Longitude"
        labelY = "Latitude"
    else:
        raise RuntimeError("Unknown CoordinateType")
except:
    self.logger.warn("Could not determine coordinate origin. Setting to (0.0, 0.0)")
    originX,originY = 0.0, 0.0

f=open(FILENAME,'r')

longitude=np.linspace(originX, originX+spacingX*NX, NX, endpoint=True)
latitude=np.linspace(originY, originY-spacingY*NY, NY, endpoint=True)
DATA=np.fromfile(FILENAME, dtype=datatype).reshape(NY, NX)
# flip data in y-direction since ER Mapper stores data bottom up
DATA=np.flipud(DATA)

x,y=np.meshgrid(longitude, latitude)
plt.figure(figsize=(6*(spacingX*NX/(spacingY*NY))+1, 6), dpi=100)
plt.pcolor(x, y, DATA)
locs,_=plt.xticks()
plt.xticks(locs, list(map(lambda x:"%g"%x, locs)))
plt.xlabel(labelX)
plt.ylabel(labelY)
plt.axis('tight')
plt.title(FILENAME)
plt.colorbar()

plt.show()
plt.savefig("ermapper_plot.png")

