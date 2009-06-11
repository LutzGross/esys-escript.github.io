
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

# You can shorten the execution time by reducing variable tend from 60 to 0.5

# Importing all the necessary modules required.
from esys.escript import *
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from numpy import identity,zeros,ones
import os

########################################################
# subroutine: cbphones
# Allows us to record the values of a PDE at various 
# specified locations in the model.
# Arguments:
#	domain  : domain of model
#	U       : Current time state displacement solution.
#	phones  : Geophone Locations
#	dim     : model dimesions
#	savepath: where to output the data files local is default
########################################################
def cbphones(domain,U,phones,dim,savepath=""):
   #find the number of geophones
   nphones = len(phones)
   u_pot = zeros([nphones,dim],float)
   
   for i in range(0,nphones):
     # define the location of the phone source 
     L=Locator(domain,numpy.array(phones[i]))
     # find potential at point source.
     temp = L.getValue(U)
     for j in range(0,dim):
       u_pot[i,j]=temp[j]

   # open file to save displacement at point source
   return u_pot