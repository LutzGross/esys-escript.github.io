########################################################
#
# Copyright (c) 2009-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009-2010 by University of Queensland
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

import numpy as np

# Extract the X and Y coordinates as two numpy arrays from escript coordiantes from .getX function
def toXYTuple(coords):
    coords = np.array(coords.toListOfTuples()) #convert to array
    coordX = coords[:,0]; coordY = coords[:,1] #X and Y components.
    return coordX,coordY
