
##############################################################################
#
# Copyright (c) 2009-2017 by The University of Queensland
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

from __future__ import division, print_function

__copyright__="""Copyright (c) 2009-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
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




