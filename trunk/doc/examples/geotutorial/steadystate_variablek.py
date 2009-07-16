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

# import tools
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
# dimensions:
L0=1.;L1=1.;
# height of k interface:
H=L1*0.75
# bottom temperature:
T_bot=100
# location, size and value of heat source
xc=[0.3,0.4]; r=0.1; Qc=3000
# two values for k
k0=1; k1=10
# create domain:
mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
x=mydomain.getX()
# set variable k
k=k0+(k1-k0)*wherePositive(x[1]-H)
# boundary temperature
T_D=T_bot/L1*(L1-x[1])
# heat source
Q=Qc*whereNegative(length(x-xc)-r)
# create PDE and set coefficients:
mypde=LinearPDE(mydomain)
mypde.setSymmetryOn()
# set PDE coefficients:
mypde.setValue(A=k*kronecker(mydomain),Y=Q, r=T_D, \
                q=whereZero(x[1])+whereZero(x[1]-L1))
# get temperature:
T=mypde.getSolution()
# save as VTK for visualisation:
saveVTK("u.vtu",T=T)

