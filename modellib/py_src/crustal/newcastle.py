"""
mining activities in modelframe 
                                                                                                                                                                                                     
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

d=True
from setups import MiningHistory
from esys.modellib.geometry import FinleyReader,VectorConstrainerOverBox
from esys.modellib.input import Sequencer
from esys.escript.modelframe import Link,Simulation, DataSource

dom=FinleyReader(debug=d)
dom.source=DataSource("./newcastle_mines.msh","gmsh")
dom.region_tag_map_source=DataSource("./vtag.xml", "ESysXML")
dom.surface_tag_map_source=DataSource("./vtag.xml", "ESysXML")

sq=Sequencer(debug=d)
sq.t=1840.
sq.t_end=2000.
sq.dt_max=100.

hist=MiningHistory(debug=d)
hist.history=DataSource("./newcastle_mining.xml")
hist.t=Link(sq,"t")

# hist.mine_locations=Link(dom,"region_tag_map")
# hist.domain=Link(dom,"domain")

fix=VectorConstrainerOverBox(debug=d)
fix.domain=Link(dom,"domain")
fix.value=0.
fix.top=False
fix.bottom=True
fix.front=False
fix.back=False
fix.left=False
fix.right=False


s=Simulation([sq, hist, fix], debug=d)
s.run()
