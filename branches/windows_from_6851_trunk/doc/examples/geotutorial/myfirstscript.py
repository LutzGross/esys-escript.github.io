from __future__ import division, print_function
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

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import sys

# get the tools we want to use
from esys.escript import *
from esys.weipa import saveVTK
try:
    from esys.dudley import Rectangle
    # some parameters
    L0=1.
    L1=1. 
    T_bot=100
    # generate n0 x n1 elements over [0,l0] x [0,l1]
    mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
    # print spatial dimension:
    print("dimension = ",mydomain.getDim())
    # get coordinates of points in domain:
    x=mydomain.getX()
    print(x) 
    # set a function 
    T_D=T_bot/L1*(L1-x[1])
    # save T_D for visualisation
    saveVTK("u.vtu",T=T_D)
except ImportError:
    print("Dudley module not available")


