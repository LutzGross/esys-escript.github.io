"""
a simple 1x1 quad`
 
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
from esys.finley import MakeDomain
p0=Point(0.,0.,0.)
p1=Point(1.,0.,0.)
p2=Point(1.,1.,0.)
p3=Point(0.,1.,0.)
l01=Line(p0,p1)
l12=Line(p1,p2)
l23=Line(p2,p3)
l30=Line(p3,p0)
c=CurveLoop(l01,l12,l23,l30)
s=PlaneSurface(c)
ps=PropertySet("The whole domain",s)
pl1=PropertySet("sides",l01,l23)
pl2=PropertySet("top and bottom",l12,l30)
d=Design(dim=2,element_size=0.2)
d.addItems(pl1,pl2)
d.addItems(ps)
d.setScriptFileName("quad.geo")
d.setMeshFileName("quad.msh")
dom=MakeDomain(d,integrationOrder=-1, reducedIntegrationOrder=-1, optimizeLabeling=True)
dom.write("quad.fly")
d.getTagMap().writeXML(open("quad_tags.xml","w")) 
# recover tagmap
#   m=TagMap()
#   m.fillFromXML(open("quad_tags.xml","r"))
