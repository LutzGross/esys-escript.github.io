
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
from esys.escript.linearPDEs import LinearPDE
import esys.dudley as dudley
from esys.weipa import saveVTK

press0=1.
lamb=1.
nu=0.3

# this sets the hook tensor:
def setHookTensor(w,l,n):
   C=Tensor4(0.,w)
   for i in range(w.getDim()):
     for j in range(w.getDim()):
          C[i,i,j,j]+=l
          C[j,i,j,i]+=n
          C[j,i,i,j]+=n

   return C
  
# generate mesh: here 10x20 mesh of order 2
domain=dudley.Rectangle(10,20,1,l0=0.5,l1=1.0)
# get handel to nodes and elements:
e=Function(domain)
fe=FunctionOnBoundary(domain)
n=ContinuousFunction(domain)
#
# set a mask msk of type vector which is one for nodes and components set be a constraint:
#
msk=whereZero(n.getX()[0])*[1.,1.]
#
#  set the normal stress components on face elements.
#  faces tagged with 21 get the normal stress [0,-press0].
#
# now the pressure is set to zero for x0 coordinates bigger then 0.1
press=whereNegative(fe.getX()[0]-0.1)*200000.*[1.,0.]
# assemble the linear system:
mypde=LinearPDE(domain)
mypde.setValue(A=setHookTensor(e,lamb,nu),y=press,q=msk,r=[0,0])
mypde.setSymmetryOn()
mypde.getSolverOptions().setVerbosityOn()
mypde.getSolverOptions().setPreconditioner(mypde.getSolverOptions().AMG)
# solve for the displacements:
u_d=mypde.getSolution()
# get the gradient and calculate the stress:
g=grad(u_d)
stress=lamb*trace(g)*kronecker(domain)+nu*(g+transpose(g))
# write the hydrostatic pressure:
saveVTK("result.vtu",displacement=u_d,pressure=trace(stress)/domain.getDim())
