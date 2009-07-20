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

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# get the tools we want to use
from esys.escript import *
from esys.finley import Rectangle
# some parameters
L0=1.
L1=1. 
T_bot=100
# generate n0 x n1 elements over [0,l0] x [0,l1]
mydomain=Rectangle(l0=L0,l1=L1,n0=20,n1=20)
# print spatial dimension:
print "dimension = ",mydomain.getDim()
# get coordinates of points in domain:
x=mydomain.getX()
print x  
# set a function 
T_D=T_bot/L1*(L1-x[1])
# save T_D for visualisation
saveVTK("u.vtu",T=T_D)

