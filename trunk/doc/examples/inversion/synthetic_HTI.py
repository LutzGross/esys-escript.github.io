from __future__ import print_function
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
from esys.downunder import Ricker, HTIWave, SimpleSEGYWriter
from math import ceil

DIM=2          # spatial dimension

ne_z=550.

absorption_zone=500*U.m

# layers from the bottom up:
layers=[absorption_zone, 100*U.m, 200 * U.m,  300*U.m ]
v_Ps=[ 3*U.km/U.sec, 3*U.km/U.sec, 1.5*U.km/U.sec, 1.5*U.km/U.sec ]
v_Ss=[1.4*U.km/U.sec, 1.4*U.km/U.sec, 0.7*U.km/U.sec, 0.7*U.km/U.sec ]
rhos=[2000*U.kg/U.m**3, 2000*U.kg/U.m**3, 2000*U.kg/U.m**3,  2000*U.kg/U.m**3]
epss=[0., 0., 0.1*0, 0. ]
gammas=[0.,0.,  0.03*0, 0.]
deltas=[0.,0., -0.1*0, 0.]        
src_dir=[0,0,1]

t_end=3.0*U.sec
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
depth=sum(layers)
#
# create array 
#
receiver_line=[2*absorption_zone + i * (rangeRcv/(numRcvPerLine-1)) for i in range(numRcvPerLine) ]
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
if DIM == 2:    
    src_locations  = [ (receiver_line[srcEW], depth) ]
    src_loc_2D=(receiver_line[srcEW], 0.)
else:
    src_locations  = [ (receiver_line[srcEW], receiver_line[srcNS], depth)]
    src_loc_2D=(receiver_line[srcEW], receiver_line[srcNS])

#
#   create sensor arrays:
#
# East-west line of receiver
rcv_locations=[]
rg=[]
mid_point=receiver_line[len(receiver_line)/2]

for ix in range(len(receiver_line)):
        if DIM == 2:
            rcv_locations.append((receiver_line[ix],  depth))
            rg.append( ( receiver_line[ix], 0.) ) 
        else:
           rcv_locations.append((receiver_line[ix], mid_point, depth))
           rg.append( ( receiver_line[ix], mid_point) ) 
# North-south line of receiver
if DIM == 3:
     for iy in range(len(receiver_line)):
            rcv_locations.append((mid_point, receiver_line[iy],  depth))
            rg.append( (  mid_point, receiver_line[iy]) ) 
#
# create domain:
#
if DIM == 2:
   domain=Rectangle(ceil(ne_z*width_x/depth), ne_z ,l0=width_x, l1=depth, 
        diracPoints=src_locations, diracTags=src_tags)
else:
   domain=Brick(ceil(ne_z*width_x/depth),ceil(ne_z*width_y/depth),ne_z,l0=width_x,l1=width_y,l2=depth, 
        diracPoints=src_locations, diracTags=src_tags)
wl=Ricker(frq)

#======================================================================
z=Function(domain).getX()[DIM-1]
z_bottom=0
v_p=0
v_s=0
delta=0
vareps=0
gamma=0
rho=0
for l in range(len(layers)):
       m=wherePositive(z-z_bottom)*whereNonPositive(z-(z_bottom+layers[l]))
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
tracer_z=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='z-displacement')
if DIM==3:
   tracer_y=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='y-displacement')

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
    print(t, loc(u[0])[len(rg)/2-4:len(rg)/2+1], wl.getValue(t))
    #if n%5 == 0 : saveSilo("tmp/u_%d.silo"%(n/5,), u=u)
    if t>0.3 and t< 0.5: 
        saveSilo("tmp/u_%d.silo"%(k,), u=u)
        k+=1
        n+=1
tracer_x.write('line_x.sgy')
tracer_z.write('line_z.sgy')
if DIM == 3: 
        tracer_y.write('line_y.sgy')

