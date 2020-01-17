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
from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript import unitsSI as U
from esys.escript.pdetools import Locator
try:
    from esys.weipa import saveSilo
    HAVE_SILO = True
except:
    HAVE_SILO = False
from esys.downunder import Ricker, SimpleSEGYWriter, HTIWave
from math import ceil
from time import time

try:
    from esys.speckley import Brick, Rectangle
    HAVE_SPECKLEY=True
except ImportError:
    HAVE_SPECKLEY=False

if HAVE_SPECKLEY and HAVE_SILO:
    DIM=2          # spatial dimension


    ORDER = 5
    ne_z= 20

    # layers from the bottom up:
    layers=[20*U.m, 180*U.m ]
    v_Ps=[i*U.km/U.sec for i in [3, 2.5]]
    v_Ss= [i*U.km/U.sec for i in [3, 2]]
    rhos=[i*U.kg/U.m**3 for i in [2.6, 2.1]]
    epss=[0, .110]
    gammas=[0, 0.035]
    deltas=[0, 0.255]
    src_dir=[0,0,1]

    t_end=0.01*U.sec #increase this end time as desired
    frq=50.*U.Hz
    sampling_interval=2*U.msec
    numRcvPerLine=101
    rangeRcv=200*U.m


    # location of source
    if DIM == 2:
        src_locations = [(0, 0)]
    else:
        src_locations = [(0, 0, 0)]

    # domain dimensions
    width_x=rangeRcv
    width_y=width_x
    depth=sum(layers)
    #
    # create array
    #
    receiver_line=[i * (rangeRcv/(numRcvPerLine-1)) for i in range(numRcvPerLine)]
    #
    #   set source location with tag "source""
    #
    src_tags=["source"]

    src_loc_2D=(0, 0)



    #
    #   create sensor arrays:
    #
    # East-west line of receivers
    rcvEW_locations=[]
    # North-south line of receivers (if 3 dimensional problem)
    rcvNS_locations=[]
    rgEW=[]
    rgNS=[]
    mid_point=receiver_line[len(receiver_line)//2]

    for ix in range(len(receiver_line)):
        rgEW.append((receiver_line[ix], 0))
        if DIM == 2:
            rcvEW_locations.append((receiver_line[ix],  0))
        else:
            rcvEW_locations.append((receiver_line[ix], 0, 0))
            rcvNS_locations.append((0, receiver_line[ix], 0))
            rgNS.append((0, receiver_line[ix]))
    # North-south line of receivers
    if DIM == 3:
         for iy in range(len(receiver_line)):
                rcv_locations.append((mid_point, receiver_line[iy],  0))
                rg.append( (  mid_point, receiver_line[iy]) )
    #
    # create domain:
    #
    if DIM == 2:
        domain = Rectangle(ORDER,
                ceil(ne_z*width_x/depth), ne_z ,l0=width_x, l1=(-depth,0),
                diracPoints=src_locations, diracTags=src_tags)
        #suppress the x-component on the x boundary
        q = whereZero(domain.getX()[0])*[1,0]
    else:
        domain=Brick(ORDER,
                ceil(ne_z*width_x/depth), ceil(ne_z*width_y/depth), ne_z,
                l0=width_x, l1=width_y, l2=(-depth,0),
                diracPoints=src_locations, diracTags=src_tags)
        q = wherePositive(
                #suppress the x-component on the x boundary
                whereZero(domain.getX()[0])*[1,0,0]
                + #logical or
                #suppress the y-component on the y boundary at the source
                whereZero(domain.getX()[1])*[0,1,0])

    # set up reciever locations
    locEW=Locator(domain,rcvEW_locations)
    tracerEW_x=SimpleSEGYWriter(receiver_group=rgEW, source=src_loc_2D,
            sampling_interval=sampling_interval,
            text='x-displacement - east-west line')
    tracerEW_z=SimpleSEGYWriter(receiver_group=rgEW, source=src_loc_2D,
            sampling_interval=sampling_interval,
            text='z-displacement - east-west line')
    if DIM==3:
        locNS=Locator(domain,rcvNS_locations)
        tracerEW_y=SimpleSEGYWriter(receiver_group=rgEW, source=src_loc_2D,
                sampling_interval=sampling_interval,
                text='x-displacement - east-west line')
        tracerNS_x=SimpleSEGYWriter(receiver_group=rgNS, source=src_loc_2D,
                sampling_interval=sampling_interval,
                text='x-displacement - north-south line')
        tracerNS_y=SimpleSEGYWriter(receiver_group=rgNS, source=src_loc_2D,
                sampling_interval=sampling_interval,
                text='y-displacement - north-south line')
        tracerNS_z=SimpleSEGYWriter(receiver_group=rgNS, source=src_loc_2D,
                sampling_interval=sampling_interval,
                text='z-displacement - north-south line')
    if not tracerEW_x.obspy_available():
        print("\nWARNING: obspy not available, SEGY files will not be written\n")
    elif getMPISizeWorld() > 1:
        print("\nWARNING: SEGY files cannot be written with multiple processes\n")


    #======================================================================
    z=ReducedFunction(domain).getX()[DIM-1]
    z_bottom=-depth
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

    wl=Ricker(frq)
    dt=min((1./5.)*min(inf(domain.getSize()/v_p), inf(domain.getSize()/v_s)), wl.getTimeScale())

    sw=HTIWave(domain, v_p, v_s, wl, src_tags[0], source_vector = src_dir,
            eps=vareps, gamma=gamma, delta=delta, rho=rho,
            absorption_zone=None, absorption_cut=1e-2, lumping=True, dt=dt)
    sw.setQ(q)

    locEW=Locator(domain, rcvEW_locations)
    if DIM == 3:
        locNS=Locator(domain, rcvNS_locations)

    mkDir('output')

    t=0.
    n=0
    k=0
    u=None
    while t < t_end:
        start = time()
        t,u = sw.update(t+sampling_interval)
        tracerEW_x.addRecord(locEW(u[0]))
        tracerEW_z.addRecord(locEW(u[DIM-1]))
        if DIM==3:
               tracerEW_y.addRecord(locEW(u[1]))
               tracerNS_x.addRecord(locNS(u[0]))
               tracerNS_y.addRecord(locNS(u[1]))
               tracerNS_z.addRecord(locNS(u[2]))
        print(t, locEW(u[DIM-1])[len(rgEW)//2-4:len(rgEW)//2+1], wl.getValue(t))
        k+=1
        if k%5 == 0:
            try:
                saveSilo("output/normalHTI_%d.silo"%(n,), v_p=v_p, u=u, cycle=k, time=t)
            except:
                print("Failed saving silo file. Was escript build without Silo support?")
            n += 1
    if k%5 != 0:
        try:
            saveSilo("output/normalHTI_%d.silo"%(n,), v_p=v_p, u=u, cycle=k, time=t)
        except:
            print("Failed saving silo file. Was escript build without Silo support?")
    if tracerEW_x.obspy_available() and getMPISizeWorld() == 1:
        tracerEW_x.write('output/lineEW_x.sgy')
        tracerEW_z.write('output/lineEW_z.sgy')
        if DIM == 3: 
            tracerEW_y.write('output/lineEW_y.sgy')
            tracerNS_x.write('output/lineNS_x.sgy')
            tracerNS_y.write('output/lineNS_y.sgy')
            tracerNS_z.write('output/lineNS_z.sgy')

else: # no speckley
    print("The Speckley module is not available")

