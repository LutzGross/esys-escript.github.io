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

from esys.pyvisi import Scene, DataCollector, ContourOnPlaneClip, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create one contour on plance clip instance.
ctropc1 = ContourOnPlaneClip(scene = s, data_collector = dc1, scalar = None, 
        viewport  = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
ctropc1.setPlaneToXY()
ctropc1.rotateY(20)
ctropc1.generateContours(8)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.elevation(-40)

s.render()

