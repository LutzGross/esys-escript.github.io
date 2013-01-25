
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
mining activities in modelframe 

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

d=True
from setups import MiningHistory, DensityChange, LinearElasticStressChange, CoulombFailureStress
from esys.modellib.geometry import FinleyReader,VectorConstrainerOverBox
from esys.modellib.input import Sequencer
from esys.escript.modelframe import Link,Simulation, DataSource
import numpy
from esys.modellib.visualization import WriteVTK

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
fix.top=[False, False, False]
fix.bottom= [True, True, True]
fix.front= [False, False, False]
fix.back=[False, False, False]
fix.left=[False, False, False]
fix.right=[False, False, False]

el=LinearElasticStressChange(debug=d)
el.domain=Link(dom,"domain")
el.tag_map=Link(dom,"tag_map")
el.density=8e3*0
el.lame_lambda=1.7e11
el.lame_mu=1.7e11
el.location_of_fixed_displacement=Link(fix,"location_of_constraint")
el.density_rate=Link(dens_dot,"density_rate")

cfs=CoulombFailureStress(debug=d)
cfs.stress=Link(el,"stress")
cfs.friction_coefficient=0.
cfs.normal=numpy.array([-1,0,1])

vis=WriteVTK()
vis.t=Link(sq)
vis.data0=Link(el,"displacement")
vis.data1=Link(cfs,"cfs")
vis.dt=10.
vis.filename="out.xml"


s=Simulation([sq, hist, dens_dot, fix, el, vis], debug=d)
s.run()
