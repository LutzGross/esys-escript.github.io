#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

__copyright__="""  Copyright (c) 2006, 2007 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

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
des.addItems(s)

dom=MakeDomain(des)
dom.write("rec.fly")
