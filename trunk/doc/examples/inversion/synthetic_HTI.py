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
from esys.downunder import Ricker, HTIWave, SimpleSEGYWriter
from math import ceil


DIM=2          # spatial dimension

v_p_top=1.5*U.km/U.sec
v_p_bottom=3*U.km/U.sec
absorption_zone=300*U.m
ne_z=500.

if DIM==3:
        # layers from the bottom up:
        layers=[500*U.m, 500 * U.m ]
        v_Ps=[ 3*U.km/U.sec, 2.*U.km/U.sec]
        v_Ss=[1.3*U.km/U.sec, 0.9*U.km/U.sec ]
        rhos=[rho=2000*U.kg/U.m**3, rho=2000*U.kg/U.m**3]
        epss=[0., 0.1]
        gammas=[0., 0.03]
        deltas=[0.,  0.1]

        depth=sum(layers)
        
        src_dir=[0,0,1]
else:
        # this is for DIM=2 case:
        v_p=2*U.km/U.sec
        v_s=0.9*U.km/U.sec
        rho=2000*U.kg/U.m**3

        vareps=0.1
        gamma=0.15
        #delta=-0.1
        delta=0.0
        delta=0.1
        depth=1*U.km
        src_dir=[1,1,0]

t_end=0.8*U.sec
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
receiver_line=[2*absorption_zone + i * (rangeRcv/(numRcvPerLine-1)) for i in xrange(numRcvPerLine) ]
#
#   set source location with tag "source""
#
src_tags=["source"]

if srcEW:
      srcNS=numRcvPerLine/2
elif srcNS:
      srcEW=numRcvPerLine/2
else:
    raise ValueError("on of the variables srcEW or srcNS must be None!")
src_loc_2D=(receiver_line[srcEW], receiver_line[srcNS])
if DIM == 2:    
    src_locations  = [ src_loc_2D ]
else:
    src_locations  = [ (receiver_line[srcEW], receiver_line[srcNS], depth)]
#
#   create sensor arrays:
#
# East-west line of receiver
rcv_locations=[]
rg=[]
mid_point=receiver_line[len(receiver_line)/2]

for ix in xrange(len(receiver_line)):
        if DIM == 2:
            rcv_locations.append((receiver_line[ix], mid_point))
            rg.append( ( receiver_line[ix], mid_point) ) 
        else:
           rcv_locations.append((receiver_line[ix], mid_point, depth))
           rg.append( ( receiver_line[ix], mid_point) ) 
# North-south line of receiver
for iy in xrange(len(receiver_line)):
        if DIM == 2:
            rcv_locations.append((mid_point, receiver_line[iy])) 
            rg.append( ( mid_point, receiver_line[iy]) ) 
        else:
            rcv_locations.append((mid_point, receiver_line[iy],  depth))
            rg.append( (  mid_point, receiver_line[iy]) ) 
#
# create domain:
#
if DIM == 2:
   domain=Rectangle(ceil(ne_z*width_x/depth),ceil(ne_z*width_y/depth),l0=width_x,l1=width_y, 
		diracPoints=src_locations, diracTags=src_tags)
else:
   domain=Brick(ceil(ne_z*width_x/depth),ceil(ne_z*width_y/depth),ne_z,l0=width_x,l1=width_y,l2=depth, 
		diracPoints=src_locations, diracTags=src_tags)
wl=Ricker(frq)

#======================================================================
if DIM == 3:
   z=Function(domain).getX()[DIM-1]
   z_bottom=0
   v_p=0
   delta=0
   vareps=0
   gamma=0
   rho=0
   for l in xrange(len(layers))
       m=wherePositive(z-z_bottom)*whereNonNegative(z-(z_bottom+layers[l]))
       v_p=v_p*(1-m)+v_Ps[l]*m
       v_s=v_s*(1-m)+v_Ss[l]*m
       rho=rho*(1-m)+rhos[l]*m
       vareps=vareps*(1-m)+epss[l]*m
       gamma=gamma*(1-m)+gammas[l]*m
       delta=delta*(1-m)+deltas[l]*m
       z_bottom+=layers[l]

sw=HTIWave(domain, v_p, v_s, wl, src_tags[0], source_vector = src_dir, eps=vareps, gamma=gamma, delta=delta, rho=rho,  \
                     absorption_zone=300*U.m, absorption_cut=1e-2, lumping=True)

loc=Locator(domain,rcv_locations)
tracer_x=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='x-displacement')
tracer_y=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='y-displacement')
if DIM==3:
   tracer_z=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='z-displacement')

t=0.
mkDir('tmp')
n=0
k=0
while t < t_end:
	t,u = sw.update(t+sampling_interval)
	tracer_x.addRecord(loc(u[0]))
	tracer_y.addRecord(loc(u[1]))
	if DIM==3:
	       tracer_z.addRecord(loc(u[2]))
	print t, loc(u[0])[len(rg)/2-4:len(rg)/2+1], wl.getValue(t)
	#if n%5 == 0 : saveSilo("tmp/u_%d.silo"%(n/5,), u=u)
	if t>0.3 and t< 0.5: 
	   saveSilo("tmp/u_%d.silo"%(k,), u=u)
	   k+=1
        n+=1
tracer_x.write('line_x.sgy')
tracer_y.write('line_y.sgy')
if DIM == 3: 
        tracer_z.write('line_z.sgy')

