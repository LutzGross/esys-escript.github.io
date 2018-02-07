
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

from esys.escript import *
from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.finley import MakeDomain


p0=Point(0.,0.)
p1=Point(1.,0.)
p2=Point(1.,1.)
p3=Point(0.,1.)

l01=Line(p0,p1)
l12=Line(p1,p2)
l23=Line(p2,p3)
l30=Line(p3,p0)

s=PlaneSurface(CurveLoop(l01,l12,l23,l30))
des=Design(dim=2, order=1, element_size = 1, keep_files=True)
des.setMeshFileName("rec.geo")
des.addItems(s)

dom=MakeDomain(des)
dom.write("rec.fly")
