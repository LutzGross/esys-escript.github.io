########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

from esys.pycad import *
import numpy as np
# To run the function it is necessary to import the modules we 
# require.
# This imports everything from the escript library
from esys.escript import * 
# numpy for array handling
import numpy as np
# pylab for matplotlib and plotting
import pylab as pl
# tools for dealing with PDEs - contains locator
from esys.escript.pdetools import *

# routine to find consecutive coordinates of a loop in pycad
def getLoopCoords(loop):
	# return all construction points of input
	temp = loop.getConstructionPoints()
	#create a numpy array for xyz components or construction points
	coords = np.zeros([len(temp),3],float)
	#place construction points in array
	for i in range(0,len(temp)):
		coords[i,:]=temp[i].getCoordinates()
	#return a numpy array
	return coords
	
# Calculate the location of quivers for a matplotlib plot
# quivshape :: [x,y] :: number of quivers in x and y direction
# lenxax :: length of model along x
# lenyax :: length of model along y
# qu :: gradient of escript solution ie grad(T)
def toQuivLocs(quivshape,lenxax,lenyax,qu):
    numquiv = quivshape[0]*quivshape[1] # total number of quivers
    dx = lenxax/quivshape[0]+1. # quiver x spacing
    dy = lenyax/quivshape[1]+1. # quiver y spacing
    qulocs=np.zeros([numquiv,2],float) # memory for quiver locations
    # fill qulocs
    for i in range(0,quivshape[0]-1):
    	for j in range(0,quivshape[1]-1):
			qulocs[i*quivshape[0]+j,:] = [dx+dx*i,dy+dy*j]
    # retreive values for quivers direction and shape from qu
    quL = Locator(qu.getFunctionSpace(),qulocs.tolist())
    qu = quL(qu) #array of dx,dy for quivers
    qu = np.array(qu) #turn into a numpy array
    qulocs = quL.getX() #returns actual locations from data
    qulocs = np.array(qulocs) #turn into a numpy array
    return qu,qulocs	

# Extract the X and Y coordinates of an array
# coords :: escript coordiantes from .getX function
def toXYTuple(coords):
    coords = np.array(coords.toListOfTuples()) #convert to Tuple
    coordX = coords[:,0]; coordY = coords[:,1] #X and Y components.
    return coordX,coordY
