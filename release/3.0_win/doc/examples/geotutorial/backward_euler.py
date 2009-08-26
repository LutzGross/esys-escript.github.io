########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
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
# end of simulation time
t_end=0.1
# time step size:
dt=0.01
# dimensions:
L0=1.;L1=1.
# location, size and value of heat source
xc=[0.3,0.4]; r=0.1; Qc=3000
# material parameter
k=1; rhocp=100; 
# bottom temperature:
T_bot=100
# generate domain:
mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
x=mydomain.getX()
# set boundray temperature:
T_D=T_bot/L1*(L1-x[1])
# set heat source:
Q=Qc*whereNegative(length(x-xc)-r)
# generate domain:
mypde=LinearPDE(mydomain)
mypde.setSymmetryOn()
# set PDE coefficients:
mypde.setValue(A=dt*k*kronecker(mydomain), D=dt*rhocp, 
                r=T_D, q=whereZero(x[1])+whereZero(x[1]-L1))
# initial temperature
T=T_D 
# step counter and time marker:
N=0; t=0
# stop when t_end is reached:
while t<t_end:
    print N,"-th time step t=",t
    # update PDE coefficient:
    mypde.setValue(Y=dt*rhocp*T+dt*Q)
    # new temperature:
    T=mypde.getSolution()
    # save as VTK for visualisation:
    saveVTK("u.%s.vtu"%N,T=T)
    # increase counter and marker:
    N+=1; t+=dt
