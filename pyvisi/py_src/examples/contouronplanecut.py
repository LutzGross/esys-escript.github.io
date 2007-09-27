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

from esys.pyvisi import Scene, DataCollector, ContourOnPlaneCut, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1000, 
        y_size = 800)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0271.vtk")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create one contour on plane cut instance.
ctropc1 = ContourOnPlaneCut(scene = s, data_collector = dc1, scalar = None, 
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
ctropc1.generateContours(8)
ctropc1.setPlaneToXY(200000)

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.elevation(-45)

s.render()

