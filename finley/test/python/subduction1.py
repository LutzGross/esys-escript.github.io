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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

# $Id$

from esys.escript import *
from esys.escript import unitsSI as U
from esys.escript.pdetools import MaskFromBoundaryTag
from esys.finley import ReadMesh
from esys.weipa import saveVTK
from esys.escript.models import StokesProblemCartesian
from math import ceil
#
#  Parameter
#
DIM=2
MESHFILE="sub.fly"
ETA=1.e22*U.Pa*U.sec
V_MAX=1.*U.cm/U.yr
ALPHA=30*U.DEG
STRIKE=10*U.DEG
DIP=30*U.DEG
N=1  # boudary layer control


g=9.81*U.m/U.sec**2
#
#  derived values
#
dom=ReadMesh(MESHFILE)
DIM=dom.getDim()
bb=boundingBox(dom)
LX=bb[0][1]-bb[0][0]
if DIM == 3: LY=bb[1][1]-bb[1][0]
DEPTH=bb[DIM-1][1]-bb[DIM-1][0]

sc=StokesProblemCartesian(dom)
x = dom.getX()
#
v=Vector(0.,Solution(dom))
mask=Vector(0.,Solution(dom))
#
#  in subduction zone:
#

if DIM==2: 
    S=numpy.array([0.,0.])
    X0=[bb[0][0],bb[1][1]]
    dd=[-cos(ALPHA),-sin(ALPHA)]
else:
    S=numpy.array([sin(STRIKE),cos(STRIKE),0.])
    X0=[bb[0][0],bb[1][0],bb[2][1]]
    dd=[-cos(ALPHA),0.,-sin(ALPHA)]
r=sqrt(length(x-X0)**2-inner(X0-x,S)**2)
v=V_MAX*r*dd
mask=MaskFromBoundaryTag(dom,"subduction")*[ 1. for i in range(DIM) ]
#
#  back of the domain
#
v=v*(1.-whereZero(x[0]-bb[0][1])*kronecker(DIM)[0])
mask+=whereZero(x[0]-bb[0][1])*kronecker(DIM)[0]
#
#  bottom of the domain
#
v=v*(1.-((bb[DIM-1][1]-x[DIM-1])/DEPTH)**N*kronecker(DIM)[DIM-1])
mask+=whereZero(x[DIM-1]-bb[DIM-1][0])*kronecker(DIM)[DIM-1]
#
#  faces of the domain:
#
if DIM==3:
    v=v*(1.-(((x[1]-bb[1][0])*(bb[1][1]-x[1])/(0.5*LY)**2)**N*kronecker(DIM)[1]))
    mask+=(whereZero(x[1]-bb[1][0])+whereZero(x[1]-bb[1][1]))*kronecker(DIM)[1]
sc.initialize(eta=ETA, fixed_u_mask= mask)
p=Scalar(0.,ReducedSolution(dom))
v,p=sc.solve(v,p, verbose=True)
saveVTK("u.vtu",velocity=v,pressure=p,m=mask)

