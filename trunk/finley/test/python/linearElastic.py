
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
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys import finley
from esys.weipa import saveVTK

pres0=-100.
lame=1.
mu=0.3
rho=1.
g=9.81

# generate mesh: here 20x20 mesh of order 1
domain=finley.Rectangle(20,20,1,l0=1.0,l1=1.0)
#
# set a mask msk of type vector which is one for nodes and components set be a constraint:
#
msk=whereZero(domain.getX()[0])*[1.,1.]
#
#  set the normal stress components on face elements.
#  faces tagged with 21 get the normal stress [0,-press0].
#
# now the pressure is set to zero for x0 coordinates equal 1. (= right face)
press=whereZero(FunctionOnBoundary(domain).getX()[0]-1.)*pres0*[1.,0.]
# assemble the linear system:
mypde=LinearPDE(domain)
k3=kronecker(domain)
k3Xk3=outer(k3,k3)

mypde.setValue(A=mu * ( swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3) ) + lame*k3Xk3, 
               Y=[0,-g*rho],
               y=press,
               q=msk,r=[0,0])
mypde.setSymmetryOn()
mypde.getSolverOptions().setVerbosityOn()
# use direct solver (default is iterative)
#mypde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
# mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)
# solve for the displacements:
u_d=mypde.getSolution()

mypde.applyOperator(u_d)
# get the gradient and calculate the stress:
g=grad(u_d)
stress=lame*trace(g)*kronecker(domain)+mu*(g+transpose(g))
# write the hydrostatic pressure:
saveVTK("result.vtu",displacement=u_d,pressure=trace(stress)/domain.getDim())
