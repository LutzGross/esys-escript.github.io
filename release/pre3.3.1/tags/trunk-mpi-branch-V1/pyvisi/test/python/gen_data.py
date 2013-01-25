# $Id:$
#
#  simple script to generate some test data for pyvisi
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
from esys.escript import *
from esys.escript.pdetools import Projector
from esys import finley
#
# 
#
DIM=2
L0=3
L1=3.
L2=1.
NE0=15
NE1=int(NE0/L0*L1)
NE2=int(NE0/L0*L2)
#
#  generate domains:
#
for DIM in [2,3]:
  if DIM == 2:
   domain=finley.Rectangle(NE0,NE1,1,l0=L0,l1=L1)
   trafo=numarray.array([[0.,1.],[-1.,0.]])
   cen=numarray.array([L0/2.,L1/2.])
  else:
   domain=finley.Brick(NE0,NE1,NE2,1,l0=L0,l1=L1,l2=L2)
   trafo=numarray.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]])
   cen=numarray.array([L0/2.,L1/2.,L2/2.])
  pp=Projector(domain)
  #
  # get function spaces:
  #
  c=ContinuousFunction(domain)
  f=Function(domain)
  b=FunctionOnBoundary(domain)
  #
  # get coordinates
  #
  c_x=c.getX()-cen
  c_r=length(c_x)
  c_t=matrix_mult(trafo,c_x/(c_r+1.e-15))
  
  f_x=f.getX()-cen
  f_r=length(f_x)
  
  b_x=b.getX()-cen
  b_r=length(b_x)
  #
  #  
  #
  saveVTK("interior_%dD.xml"%DIM,temperature=sin(2*c_r),temperature_cell=sin(2*f_r), \
                                 velocity=matrix_mult(trafo,c_x/(c_r+1.e-15)), velocity_cell=matrix_mult(trafo,f_x/f_r), \
                                 stress=pp(grad(c_t)),                                       \
                                 stress_cell=grad(c_t))
  saveVTK("boundary_%dD.xml"%DIM,temperature=sin(b_r),velocity=matrix_mult(trafo,b_x/(b_r+1.e-15)), \
                                 stress=grad(c_t,b))
