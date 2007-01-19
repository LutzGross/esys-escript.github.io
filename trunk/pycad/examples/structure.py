"""
seismic wave propagation domain
 
@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""


__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from esys.pycad import *
from esys.pycad.gmsh import Design

l=100000.           # width and length m (without obsorber)
h=30000.            # width and length m (without obsorber)
d_absorber=l*0.10   # thickness of absorbing layer

b=Brick(Point(0.,0.,-h),Point(l,l,0.))
p1=Point(l/5,l/5,-2*h/3)
p2=p1+[l/5,l/5,0.]
p3=p2+[0,0.,h/5]
p4=p3+[-l/5,-l/5,0.]
l1=Line(p1,p2)
l2=Line(p2,p3)
l3=Line(p3,p4)
l4=Line(p4,p1)
c1=CurveLoop(l1,l2,l3,l4)
c1.setLocalScale(0.1)
q1=Point(2*l/3,2*l/3,-2*h/3,local_scale=0.3)
q2=q1+[-l/4,l/4,0.]
q3=q2+[0,0.,h/3]
q4=q3+[l/4,-l/4,0.]
m1=Line(q1,q2)
m2=Line(q2,q3)
m3=Line(q3,q4)
m4=Line(q4,q1)
c2=CurveLoop(m1,m2,m3,m4)
c2.setLocalScale(0.1)
dsgn=Design(element_size=h/3)
dsgn.addItems(Volume(SurfaceLoop(*tuple(b.getSurfaces()+[PlaneSurface(c1),PlaneSurface(c2)]))))
dsgn.setScriptFileName("test.geo")
dsgn.setMeshFileName("test.msh")
print dsgn.getMeshHandler()
