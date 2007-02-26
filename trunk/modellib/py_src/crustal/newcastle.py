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
from setups import MiningHistory, DensityChange, LinearElasticStressChange
from esys.modellib.geometry import FinleyReader,VectorConstrainerOverBox
from esys.modellib.input import Sequencer
from esys.escript.modelframe import Link,Simulation, DataSource

dom=FinleyReader(debug=d)
dom.source=DataSource("./newcastle_mines.msh","gmsh")
dom.tag_map_source=DataSource("./tags.xml", "ESysXML")

sq=Sequencer(debug=d)
sq.t=1840.
sq.t_end=2000.
sq.dt_max=100.

hist=MiningHistory(debug=d)
hist.history=DataSource("./newcastle_mining.xml")
hist.t=Link(sq,"t")

dens_dot=DensityChange(debug=d)
dens_dot.domain=Link(dom,"domain")
dens_dot.tag_map=Link(dom,"tag_map")
dens_dot.mass_rate=Link(hist,"mass_changes")

fix=VectorConstrainerOverBox(debug=d)
fix.domain=Link(dom,"domain")
fix.value=0.
fix.top=False
fix.bottom=True
fix.front=False
fix.back=False
fix.left=False
fix.right=False
# dens_dot.density_rate=Link(hist,"mass_changes")

el=LinearElasticStressChange(debug=d)
el.domain=Link(dom,"domain")
el.tag_map=Link(dom,"tag_map")
el.density=1.
el.lame_lambda=2.
el.lame_mu=1.
el.location_of_fixed_displacement=Link(fix,"location_of_constraint")
el.density_rate=Link(dens_dot,"density_rate")

# hist.domain=Link(dom,"domain")


s=Simulation([sq, hist, dens_dot, fix, el], debug=d)
s.run()
