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
