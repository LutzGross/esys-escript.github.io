
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

"""
some basic shapes.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from .primitives import *

def Brick(start,end):
    """
    Creates a brick with given start and end point.
    """
    dx=end.getCoordinates()-start.getCoordinates()
    p000=start+[   0.,0.,0.]
    p100=start+[dx[0],0.,0.]
    p010=start+[0.,dx[1],0.]
    p110=start+[dx[0],dx[1],0.]
    p001=start+[0.,0.,dx[2]]
    p101=start+[dx[0],0.,dx[2]]
    p011=start+[0.,dx[1],dx[2]]
    p111=start+[dx[0],dx[1],dx[2]]
    l10=Line(p000,p100)
    l20=Line(p100,p110)
    l30=Line(p110,p010)
    l40=Line(p010,p000)
    l11=Line(p000,p001)
    l21=Line(p100,p101)
    l31=Line(p110,p111)
    l41=Line(p010,p011)
    l12=Line(p001,p101)
    l22=Line(p101,p111)
    l32=Line(p111,p011)
    l42=Line(p011,p001)
    bottom=PlaneSurface(CurveLoop(-l10,-l40,-l30,-l20))
    top=PlaneSurface(CurveLoop(l12,l22,l32,l42))
    front=PlaneSurface(CurveLoop(-l11,l10,l21,-l12))
    back=PlaneSurface(CurveLoop(l30,l41,-l32,-l31))
    left=PlaneSurface(CurveLoop(l11,-l42,-l41,l40))
    right=PlaneSurface(CurveLoop(-l21,l20,l31,-l22))
    return SurfaceLoop(bottom,top,front,back,left,right)

