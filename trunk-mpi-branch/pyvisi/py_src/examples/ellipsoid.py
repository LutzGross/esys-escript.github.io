#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

from esys.pyvisi import Scene, DataCollector, Ellipsoid, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# Create one ellipsoid instance.
e1 = Ellipsoid(scene = s, data_collector = dc1, tensor = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
e1.setScaleFactor(scale_factor = 0.1)
e1.setPhiResolution(10)
e1.setThetaResolution(10)
e1.setDimension(2,2,2)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()
c1.elevation(-20)

s.render()

