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

# $Id:$

from esys.escript import *
from esys.weipa import saveVTK
from esys.escript.models import StokesProblemCartesian

try:
    from esys.finley import Rectangle
    HAVE_FINLEY = True
except ImportError as e:
    print("Finley module required but not available")
    HAVE_FINLEY = False

if HAVE_FINLEY:
    NE=25
    dom = Rectangle(NE,NE,order=-1, useElementsOnFace=0)  # use macro elements for pressure
    x = dom.getX()
    sc=StokesProblemCartesian(dom)
    mask= (whereZero(x[0])*[1.,0]+whereZero(x[0]-1))*[1.,0] + (whereZero(x[1])*[0.,1.]+whereZero(x[1]-1))*[1.,1] 
    sc.initialize(eta=.01, fixed_u_mask= mask)
    v=Vector(0.,Solution(dom))
    v[0]+=whereZero(x[1]-1.)
    p=Scalar(0.,ReducedSolution(dom))
    v,p=sc.solve(v,p, verbose=True)
    saveVTK("u.vtu",velocity=v,pressure=p)


