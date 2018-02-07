
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

from __future__ import division, print_function

import pylab as pl
import numpy as np

s=0
e=20
x = np.linspace(s,e,1000)

y = np.exp(-1.0*(0.005*(s-x))**2)

pl.plot(x,y)
pl.title("cerjen damping")



pl.savefig("cerjen.png")


x2= np.linspace(0,0.1,333)

y2=np.exp(-50.*x2)*np.sin(40*3.14157*x2)

pl.clf()
pl.plot(x2,y2)
pl.savefig("source.png")
