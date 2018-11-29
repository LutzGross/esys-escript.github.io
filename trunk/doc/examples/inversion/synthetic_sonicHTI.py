from __future__ import division
from __future__ import print_function
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################
from __future__ import print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript import unitsSI as U
from esys.escript.pdetools import Locator
from esys.weipa import saveSilo
from esys.downunder import Ricker, SonicHTIWave, SimpleSEGYWriter
from math import ceil
import time, os

try:
    from esys.ripley import Brick, Rectangle
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False

if HAVE_RIPLEY:
    DIM=2          # spatial dimension

    # layers from the bottom up:
    layers = [ 1*U.km     , 1*U.km  ,700*U.m, 500*U.m, 800*U.m ]
    v_Ps= [ 3.8 * U.km/U.sec , 3. * U.km/U.sec, 2.5*U.km/U.sec, 1.9*U.km/U.sec, 1.5*U.km/U.sec]
    epss =[   0., 0.24, 0, 0.1, 0]
    deltas=[  0.,  0.1, 0.,0.03,0 ]
    azmths=[  0.,0.,0,  0, 0.]

    dt=0.5*U.msec

    ne_z=40

    dt=0.5*U.msec

    t_end=0.008*U.sec #only this low for testing purposes
    frq=15.*U.Hz
    tcenter=None
    sampling_interval=4*U.msec
    numRcvPerLine=101
    rangeRcv=4.*U.km
    src_dir=[0,1]
    absorption_zone=1000*U.m

    # location of source in crossing array lines with in 0..numRcvInLine one needs to be None
    srcEW=numRcvPerLine//2
    srcNS=None
    # dommain dimension
    width_x=rangeRcv + 2*absorption_zone
    width_y=width_x
    depth=sum(layers)
    ne_x=int(ceil(ne_z*width_x/depth))
    #
    # create array 
    #
    receiver_line=[  absorption_zone  + i * (rangeRcv//(numRcvPerLine-1) ) for i in range(numRcvPerLine) ]
    #
    #   set source location with tag "source""
    #
    src_tags=["source"]

    if srcEW:
          srcNS=numRcvPerLine//2
    elif srcNS:
          srcEW=numRcvPerLine//2
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
    mid_point=receiver_line[len(receiver_line)//2]

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
       domain=Rectangle(ne_x, ne_z ,l0=width_x, l1=depth, 
            diracPoints=src_locations, diracTags=src_tags)
    else:
       domain=Brick(ne_x,ne_x,ne_z,l0=width_x,l1=width_y,l2=depth, 
            diracPoints=src_locations, diracTags=src_tags)
    wl=Ricker(frq, tcenter)

    #======================================================================
    z=Function(domain).getX()[DIM-1]
    z_bottom=0
    v_p=0
    delta=0
    vareps=0
    azmth=0
    rho=0
    for l in range(len(layers)):
           m=wherePositive(z-z_bottom)*whereNonPositive(z-(z_bottom+layers[l]))
           v_p=v_p*(1-m)+v_Ps[l]*m
           vareps=vareps*(1-m)+epss[l]*m
           azmth=azmth*(1-m)+azmths[l]*m
           delta=delta*(1-m)+deltas[l]*m
           z_bottom+=layers[l]

    sw=SonicHTIWave(domain, v_p, wl, src_tags[0], dt=dt, source_vector = src_dir, eps=vareps, delta=delta, azimuth=azmth,  \
                         absorption_zone=absorption_zone, absorption_cut=1e-2, lumping=False)

    #
    #  print some info:
    #
    print("ne_x = ", ne_x)
    print("ne_z = ", ne_z)
    print("h_x = ", width_x/ne_x)
    print("h_z = ", depth/ne_z)
    print("dt = ", sw.getTimeStepSize()*1000, "msec")
    print("width_x = ", width_x)
    print("depth = ", depth)
    print("number receivers = ", numRcvPerLine)
    print("receiver spacing = ", receiver_line[1]-receiver_line[0])
    print("sampling time = ", sampling_interval*1000,"msec")
    print("source @ ", src_locations[0])
    #
    loc=Locator(domain,rcv_locations)
    tracerP=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='P')
    tracerQ=SimpleSEGYWriter(receiver_group=rg, source=src_loc_2D, sampling_interval=sampling_interval, text='Q')

    if not tracerP.obspy_available():
        print("\nWARNING: obspy not available, SEGY files will not be written\n")
    elif getMPISizeWorld() > 1:
        print("\nWARNING: SEGY files cannot be written with multiple processes\n")

    t=0.
    OUT_DIR="out%sm%smus"%(int(width_x/ne_x),int(sw.getTimeStepSize()*1000000))
    mkDir(OUT_DIR)
    n=0
    k=0
    timer1=time.time()
    while t < t_end:
        t,u = sw.update(t+sampling_interval)
        Plog=loc(u[1])
        Qlog=loc(u[0])
        tracerP.addRecord(Plog)
        tracerQ.addRecord(Qlog)
        print(t, wl.getValue(t)," :", Plog[0], Plog[srcEW], Plog[-1])
    timer1=time.time()-timer1
    print("time= %e sec; %s sec per step"%(timer1,timer1/max(sw.n,1)))

    if tracerP.obspy_available() and getMPISizeWorld() == 1:
        tracerP.write(os.path.join(OUT_DIR,'lineP.sgy'))
        tracerQ.write(os.path.join(OUT_DIR,'lineQ.sgy'))

else: # no ripley
    print("The Ripley module is not available")

