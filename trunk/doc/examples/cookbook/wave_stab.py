
##############################################################################
#
# Copyright (c) 2009-2012 by University of Queensland
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

__copyright__="""Copyright (c) 2009-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# Antony Hallam
# Acoustic Wave Equation Simulation

# Importing all the necessary modules required.
import matplotlib
matplotlib.use('agg') #It's just here for automated testing

import numpy as np
import pylab as pl

#Geometric and material property related variables.
mx = 1000. # model lenght

ndx = np.arange(0.,20.,.1)
mtim= np.zeros(len(ndx),'float')
nvel= np.arange(500.,5000.,500.)

for vel in nvel:
	mtim=ndx/vel
	pl.plot(ndx,mtim,label='%d m/s'%vel)
	
pl.title('Maximum time steps calculations by velocity')
pl.xlabel('Minimum grid spacing (m)')
pl.ylabel('Maximum stable time step (s)')
pl.legend()
pl.show()




