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

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
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
from esys.downunder import Ricker, TTIWave, SimpleSEGYWriter
from math import ceil
import time
try:
    from esys.speckley import Rectangle
    HAVE_SPECKLEY=True
except ImportError:
    HAVE_SPECKLEY=False

if HAVE_SPECKLEY and HAVE_SILO:
    # these are the layers from the top down
    layers = [     400*U.m         ,    100*U.m  ,        1.*U.km,         ]
    v_P=     [    2.86* U.km/U.sec ,    1.5 * U.km/U.sec, 2.86 * U.km/U.sec     ]
    v_S=     [     1.79 * U.km/U.sec ,    0.7* U.km/U.sec, 1.8*U.km/U.sec     ]
    eps =    [    0.               ,    0.5,               0.1             ]
    delta=   [    0.               ,    0.5            ,    0.              ]
    tilt=    [    0.               ,    0.             ,    0.              ]
    rho=     [ 2000 * U.kg/U.m**3  , 2000 * U.kg/U.m**3, 2000 * U.kg/U.m**3 ]
    #
    #   other input:
    #
    t_end=0.008*U.sec                   # only this low for testing purposes
    frq=10.*U.Hz                        # dominant frequnce in the Ricker (maximum frequence ~ 2 * frq)
    sampling_interval=4*U.msec          # sampling interval
    ne_z=None                           # number of elements in vertical direction, if none it is guessed 
    n_out = 5                         # a silo file is created every n_out's sample
    absorption_zone=100*U.m             # absorbtion zone to be added in horizontal direction to the area covered by receiver line 
                                        # and subtracted from the lowest layer.
    # defines the receiver line 
    rangeRcv=800*U.m                    # width of the receiver line
    numRcvPerLine=101                   # total number of receivers
    src_id=numRcvPerLine//2             # location of source in crossing array lines with in 0..numRcvInLine 
    lumping = True
    src_dir=[0,1]

    # domain dimension
    width_x=rangeRcv + 4*absorption_zone
    depth=sum(layers)
    if ne_z is None:
        ne_z=int(ceil(depth*(2*frq)/min(v_P)))
        if getMPISizeWorld() > 10:
            ne_z = 2*ne_z-1
    ne_x=int(ceil(ne_z*width_x/depth))
    #
    # create receiver array 
    #
    receiver_line=[2*absorption_zone + i * (rangeRcv//(numRcvPerLine-1)) for i in range(numRcvPerLine) ]
    #
    #   set source location with tag "source""
    #
    src_tags=["source"]
    src_locations = [ (receiver_line[src_id], depth)]
    srcloc=(receiver_line[src_id], 0.)
    #
    #   output
    #
    print("%s"%(time.asctime(),))
    print("ne_x = %s"%(ne_x,))
    print("ne_z = %s"%(ne_z,))
    print("width = %s m"%(width_x,))
    print("depth = %s m"%(depth, ))
    print("absorption_zone = %s m"%(absorption_zone, ))
    print("sampling interval = %s ms"%(sampling_interval/U.msec,))
    print("t_end = %s sec"%(t_end,))
    print("ricker dominant freqency = %s Hz"%(frq,))
    print("length of receiver line = %s ms"%(rangeRcv,))
    print("number of receivers = %s"%(numRcvPerLine,))
    print("first receiver location = %s m"%(receiver_line[0],))
    print("last receiver location = %s m"%(receiver_line[-1],))
    print("source location = %s m"%(src_locations[0][0],))
    print("source orientation = %s"%(src_dir,))
    print("matrix lumping = %s"%(lumping,))
    print("Layer\tV_p\tV_s\teps\tdelta\ttilt\trho")
    for i in range(len(layers)):
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s"%( layers[i], v_P[i], v_S[i], eps[i], delta[i], tilt[i], rho[i]))
    #
    # create domain:
    #
    order = 5
    domain=Rectangle(order, ne_x,ne_z, l0=width_x, l1=depth, 
                diracPoints=src_locations, diracTags=src_tags, d0=getMPISizeWorld())
    #
    # create the wavelet:
    #
    wl=Ricker(frq)
    #
    #======================================================================
    #
    #  set 
    #
    z=ReducedFunction(domain).getX()[1]
    z_top=0
    V_P=0
    V_S=0
    Delta=0
    Eps=0
    Tilt=0
    Rho=0
    z_top=depth

    for l in range(len(layers)):
        m=whereNonPositive(z-z_top)*wherePositive(z-(z_top-layers[l]))
        V_P = V_P     * (1-m)  + v_P[l]  * m
        V_S = V_S     * (1-m)  + v_S[l]  * m
        Delta = Delta * (1-m)  + delta[l]* m
        Eps = Eps     * (1-m)  + eps[l]  * m
        Tilt = Tilt   * (1-m)  + tilt[l] * m
        Rho = Rho     * (1-m)  + rho[l]  * m
        z_top-=layers[l]

    sw=TTIWave(domain, V_P, V_S, wl, src_tags[0], source_vector = src_dir,
                    eps=Eps, delta=Delta, rho=Rho, theta=Tilt,
                    absorption_zone=absorption_zone, absorption_cut=1e-2, lumping=lumping)

    srclog=Locator(domain, [ (r , depth) for r in receiver_line ] )
    grploc=[ (x[0], 0.) for x in srclog.getX() ]

    tracer_x=SimpleSEGYWriter(receiver_group=grploc, source=srcloc, sampling_interval=sampling_interval, text='x-displacement')
    tracer_z=SimpleSEGYWriter(receiver_group=grploc, source=srcloc, sampling_interval=sampling_interval, text='z-displacement')

    if not tracer_x.obspy_available():
        print("\nWARNING: obspy not available, SEGY files will not be written\n")
    elif getMPISizeWorld() > 1:
        print("\nWARNING: SEGY files cannot be written with multiple processes\n")

    t=0.
    mkDir('output')
    n=0
    k_out=0
    print("calculation starts @ %s"%(time.asctime(),))
    while t < t_end:
            t,u = sw.update(t+sampling_interval)
            tracer_x.addRecord(srclog(u[0]))
            tracer_z.addRecord(srclog(u[1]))
            print("t=%s, src=%s: \t %s \t %s \t %s"%(t, wl.getValue(t),srclog(u[1])[0], srclog(u[1])[src_id], srclog(u[1])[-1]))
            if not n_out is None and n%n_out == 0:
                print("time step %s written to file %s"%(n_out, "output/u_%d.silo"%(k_out,)))
                try:
                    saveSilo("output/u_%d.silo"%(k_out,), u=u)
                except:
                    print("Failed saving silo file. Was escript build without Silo support?")
                k_out+=1
            n+=1
    if tracer_x.obspy_available() and getMPISizeWorld() == 1:
        tracer_x.write('output/lineX.sgy')
        tracer_z.write('output/lineZ.sgy')
    print("calculation completed @ %s"%(time.asctime(),))

else: # no speckley
    print("The Speckley module is not available")

