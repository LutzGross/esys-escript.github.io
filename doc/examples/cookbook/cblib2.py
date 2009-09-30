########################################################
#
# Copyright (c) 2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Additional routines using matplotlib for cookbook examples.
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""

#Please ensure that if you need to use the agg back-end that you have chosen it before 
#calling functions in this file


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

def toRegGrid(grid,domain,newx,newy,width,depth):
	import pylab as pl
	oldspacecoords=domain.getX()
	gridT = grid.toListOfTuples(scalarastuple=False)
	cxy = Data(oldspacecoords, grid.getFunctionSpace())
	cX, cY = toXYTuple(cxy)
	xi = np.linspace(0.0,width,newx)	
	yi = np.linspace(depth,0.0,newy)
	# grid the data.
	zi = pl.matplotlib.mlab.griddata(cX,cY,gridT,xi,yi)
	return xi,yi,zi
	
def gradtoRegGrid(grid,domain,newx,newy,width,depth,dim):
	import pylab as pl
	cxy = grid.getFunctionSpace().getX()
	gridT = grid.toListOfTuples()#(scalarastuple=False)
	#cxy = Data(oldspacecoords, grid.getFunctionSpace())
	cX, cY = toXYTuple(cxy)
	xi = np.linspace(0.0,width,newx)	
	yi = np.linspace(depth,0.0,newy)
	
	gridT = np.array(gridT)
	gridT = gridT[:,dim]
	
    # grid the data.
	
	zi = pl.matplotlib.mlab.griddata(cX,cY,gridT,xi,yi)
	return xi,yi,zi
