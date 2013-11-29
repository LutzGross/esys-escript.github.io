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


DIM=2           # spatial dimension

width_x=2*U.km  # area of the domain
width_y=2*U.km
depth=1*U.km    # depth 

v_p_top=1.5*U.km/U.sec
v_p_bottom=3*U.km/U.sec
reflector_at=0.5*depth

ne_z=400

t_end=0.6*U.sec

frq=20*U.Hz
sampling_interval=4*U.msec


#
# create array 
#
grid=[ width_x/2 -300*U.m, width_x/2 -200*U.m, width_x/2 -100*U.m, width_x/2, width_x/2 + 100*U.m, width_x/2 + 200*U.m, width_x/2 + 300*U.m ]
array_points=[]
array_tags=[]
rg=[]
if DIM == 2: 
   src_tag='phone_3'
   src=grid[len(grid)/2]
   for ix in xrange(len(grid)):
        array_points.append([grid[ix], depth])
        array_tags.append('phone_%s'%(ix,))
        rg.append(grid[ix]) 
else:
   src_tag='phone_33'
   src=( grid[len(grid)/2], grid[len(grid)/2])
   for ix in xrange(len(grid)):
      for iy in xrange(len(grid)):
           array_points.append([grid[ix], grid[iy], depth])
           array_tags.append('phone_%s%s'%(ix,iy))
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
m=whereNegative(domain.getX()[DIM-1]-reflector_at)
v_p=v_p_bottom*m+v_p_top*(1-m)

sw=SonicWave(domain, v_p, source_tag=src_tag, wavelet=wl, absorption_zone=300*U.m)

loc=Locator(domain,array_points)

tracer=SimpleSEGYWriter(receiver_group=rg, source=src, sampling_interval=sampling_interval)

t=0.
mkDir('tmp')
n=0
while t < t_end:
	t,p = sw.update(t+sampling_interval)
	rec=loc(p)
	tracer.addRecord(rec)
	print t, rec[:4], wl.getValue(t)
	saveSilo("tmp/u_%d.silo"%(n,), p=p)
	n+=1
tracer.write('line.sgy')

