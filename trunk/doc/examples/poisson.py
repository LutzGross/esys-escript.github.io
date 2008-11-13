
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

from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.finley import Rectangle,Brick
import time
from optparse import OptionParser

import sys

#sys.stderr=None
parser = OptionParser(usage="%prog [-e [NE]]")
parser.add_option("-e", "--elements", dest="NE", help="number of elements in one direction.",metavar="NE", default=30)
(options, args) = parser.parse_args()
NE=int(options.NE)

# generate domain:
mydomain = Rectangle(l0=1.,l1=1.,n0=NE, n1=NE)
#mydomain = Brick(l0=1.,l1=1., l2=1.,n0=10, n1=10, n2=10)
# define characteristic function of Gamma^D
x = mydomain.getX()
gammaD = whereZero(x[0])+whereZero(x[1])
# define PDE and get its solution u
mypde = Poisson(domain=mydomain)
mypde.setValue(f=1,q=gammaD)

#mypde.setSolverMethod(mypde.PCG,mypde.JACOBI)
#t1 = time.time()
#u=mypde.getSolution(verbose=True)
#elapsed = time.time() - t1
#print "ESCRIPT Job completed","\t",NE,"\tin\t",elapsed,"\t(seconds)\t",elapsed/3600.0, " (hours)"
#print u

mypde.setSolverMethod(mypde.PCG,mypde.AMG)
t1 = time.time()
u=mypde.getSolution(verbose=True)
elapsed = time.time() - t1
print "AMG Job completed","\t",NE,"\tin\t",elapsed,"\t(seconds)\t",elapsed/3600.0, " (hours)"
print u

mypde.setSolverMethod(mypde.PCG,mypde.JACOBI)
t1 = time.time()
u=mypde.getSolution(verbose=True)
elapsed = time.time() - t1
print "ESCRIPT Job completed","\t",NE,"\tin\t",elapsed,"\t(seconds)\t",elapsed/3600.0, " (hours)"
print u
#u = mypde.getSolution()
# write u to an external file
#saveVTK("u.xml",sol=u)
