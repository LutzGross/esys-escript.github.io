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
__url__="http://www.uq.edu.au/esscc/escript-finley"

# import tools
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
# set dimensions
L0=1.;L1=1. 
# bottom temperature:
T_bot=100
# location, size and value of heat source
xc=[0.3,0.4]; r=0.1; Qc=3000
# create domain
mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
x=mydomain.getX()
k=1
# temperature for boundary condition
T_D=T_bot/L1*(L1-x[1])
# heat source
Q=Qc*whereNegative(length(x-xc)-r)
# create PDE:
mypde=LinearPDE(mydomain)
mypde.setSymmetryOn()
# set coefficients:
mypde.setValue(A=k*kronecker(mydomain),Y=Q, r=T_D, \
                q=whereZero(x[1])+whereZero(x[1]-L1))
# get temperature:
T=mypde.getSolution()
# write to file:
saveVTK("u.xml",T=T)
