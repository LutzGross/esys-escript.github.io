##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################
from __future__ import print_function

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
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
absorption_zone=300*U.m
ne_z=400

reflector_at=0.5*depth


t_end=1.*U.sec
frq=20.*U.Hz
sampling_interval=4*U.msec
numRcvPerLine=101
rangeRcv=800*U.m

# location of source in crossing array lines with in 0..numRcvInLine one needs to be None
srcEW=numRcvPerLine/2
srcNS=None

# dommain dimension
width_x=rangeRcv + 4*absorption_zone
width_y=width_x
#
# create array 
#
receiver_line=[2*absorption_zone + i * (rangeRcv/(numRcvPerLine-1)) for i in range(numRcvPerLine) ]
#
#   set source location with tag "source""
#
src_tags=["source"]
if DIM == 2:
   src_locations = [ (receiver_line[srcEW], depth)]
   src_loc_2D=(receiver_line[srcEW], 0.)
else:
   if srcEW:
      srcNS=numRcvPerLine/2
   elif srcNS:
      srcNS=numRcvPerLine/2
   else:
       raise ValueError("on of the variables srcEW or srcNS must be None!")
   src_locations  = [ (receiver_line[srcEW], receiver_line[srcNS], depth)]
   src_loc_2D=(receiver_line[srcEW], receiver_line[srcNS])
#
#   create sensor arrays:
#
# East-west line of receiver
rcvEW_locations=[]
rgEW=[]
mid_point=receiver_line[len(receiver_line)/2]

for ix in range(len(receiver_line)):
        if DIM == 2:
            rcvEW_locations.append((receiver_line[ix], depth))
            rgEW.append( ( receiver_line[ix], 0.) ) 
        else:
           rcvEW_locations.append((receiver_line[ix], mid_point, depth))
           rgEW.append( ( receiver_line[ix], mid_point) ) 
# North-south line of receiver
if DIM == 3:
   rcvNS_locations=[]
   rgNS=[]

   for iy in range(len(receiver_line)):
       rcvNS_locations.append((mid_point, receiver_line[iy],  depth))
       rgNS.append( (  mid_point, receiver_line[iy]) ) 
#
# create domain:
#
if DIM == 2:
   domain=Rectangle(ceil(ne_z*width_x/depth),ne_z,l0=width_x,l1=depth, 
        diracPoints=src_locations, diracTags=src_tags)
else:
   domain=Brick(ceil(ne_z*width_x/depth),ceil(ne_z*width_y/depth),ne_z,l0=width_x,l1=width_y,l2=depth, 
        diracPoints=src_locations, diracTags=src_tags)
wl=Ricker(frq)
m=whereNegative(Function(domain).getX()[DIM-1]-reflector_at)
v_p=v_p_bottom*m+v_p_top*(1-m)

sw=SonicWave(domain, v_p, source_tag=src_tags[0], wavelet=wl, absorption_zone=absorption_zone, lumping=True)

locEW=Locator(domain,rcvEW_locations)
tracerEW=SimpleSEGYWriter(receiver_group=rgEW, source=src_loc_2D, sampling_interval=sampling_interval)
if DIM==3:
   locNS=Locator(domain,rcvNS_locations)
   tracerNS=SimpleSEGYWriter(receiver_group=rgNS, source=src_loc_2D, sampling_interval=sampling_interval)


t=0.
mkDir('tmp')
n=0
while t < t_end:
    t,p = sw.update(t+sampling_interval)
    tracerEW.addRecord(locEW(p))
    if DIM==3: tracerNS.addRecord(locNS(p))
    print(t, locEW(p)[:4], wl.getValue(t))
    if n%5 == 0 : saveSilo("tmp/u_%d.silo"%(n/5,), p=p)
    n+=1
tracerEW.write('lineEW.sgy')
if DIM == 3: tracerNS.write('lineNS.sgy')
