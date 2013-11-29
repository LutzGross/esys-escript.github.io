##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript import unitsSI as U
from esys.escript.pdetools import Locator
from esys.finley import Brick, Rectangle
from esys.weipa import saveSilo
from esys.downunder import Ricker, SonicWave, SimpleSEGYWriter
from math import ceil


DIM=2          # spatial dimension

depth=1*U.km    # depth 

v_p_top=1.5*U.km/U.sec
v_p_bottom=3*U.km/U.sec
reflector_at=0.5*depth

ne_z=400

t_end=1.0*U.sec

frq=20*U.Hz
sampling_interval=4*U.msec
if DIM == 2:
   rcv_num=100
else:
   rcv_num=30

rcv_range=600*U.m

absorption_zone=300*U.m

width_x=rcv_range + 4*absorption_zone
width_y=width_x
#
# create array 
#
grid=[2*absorption_zone + i * (rcv_range/rcv_num) for i in xrange(rcv_num) ]

array_points=[]
array_tags=[]
rg=[]
if DIM == 2:
   src_id=rcv_num/2
   for ix in xrange(len(grid)):
        array_points.append((grid[ix], depth))
        array_tags.append('phone_%3.3d'%(ix,))
        rg.append( ( grid[ix], 0.) ) 
else:
   src_id= (rcv_num/2) * rcv_num + (rcv_num/2) 
   for ix in xrange(len(grid)):
      for iy in xrange(len(grid)):
           array_points.append((grid[ix], grid[iy], depth))
           array_tags.append('phone_%3.3d-%3.3d'%(ix,iy))
           rg.append( (grid[ix], grid[iy]) )
#
# create domain:
#
if DIM == 2:
   domain=Rectangle(ceil(ne_z*width_x/depth),ne_z,l0=width_x,l1=depth, 
		diracPoints=array_points, diracTags=array_tags)
else:
   domain=Brick(ceil(ne_z*width_x/depth),ceil(ne_z*width_y/depth),ne_z,l0=width_x,l1=width_y,l2=depth, 
		diracPoints=array_points, diracTags=array_tags)
wl=Ricker(frq)
m=whereNegative(Function(domain).getX()[DIM-1]-reflector_at)
v_p=v_p_bottom*m+v_p_top*(1-m)

sw=SonicWave(domain, v_p, source_tag=array_tags[src_id], wavelet=wl, absorption_zone=absorption_zone, lumping=True)

loc=Locator(domain,array_points)

print src_id

tracer=SimpleSEGYWriter(receiver_group=rg, source=array_points[src_id], sampling_interval=sampling_interval)

t=0.
mkDir('tmp')
n=0
while t < t_end:
	t,p = sw.update(t+sampling_interval)
	rec=loc(p)
	tracer.addRecord(rec)
	print t, rec[:4], wl.getValue(t)
	if n%5 == 0 : saveSilo("tmp/u_%d.silo"%(n/5,), p=p)
	n+=1
tracer.write('line.sgy')

