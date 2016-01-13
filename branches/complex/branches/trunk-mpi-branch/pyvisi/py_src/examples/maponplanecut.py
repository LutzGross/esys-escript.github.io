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

from esys.pyvisi import Scene, DataCollector, MapOnPlaneCut, Camera
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 1152, 
        y_size = 864)

# Create two data collector instances.
dc1 = DataCollector(source = Source.XML)
dc2 = DataCollector(source = Source.XML)

dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")
dc2.setFileName(file_name = 
        "/home/jongui/data/laurent/subduction/source/function.0001.vtk")

# NOTE: There is a difference between performing rotation then followed by 
# translation and performing translation then followed by rotation.

# Create a map on plane cut instance for the first viewport.
mopc1 = MapOnPlaneCut(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST)
mopc1.setPlaneToYZ(1.5)

c1 = Camera(scene = s, data_collector = dc2, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

# Create three map on plane cut instances for the second viewport.
mopc2_1 = MapOnPlaneCut(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_WEST)
mopc2_1.setPlaneToYZ(1.5)

mopc2_2 = MapOnPlaneCut(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_WEST)
mopc2_2.setPlaneToXZ(1.5)

mopc2_3 = MapOnPlaneCut(scene = s, data_collector = dc1, 
        viewport = Viewport.NORTH_WEST)
mopc2_3.setPlaneToXY()

c2 = Camera(scene = s, data_collector = dc2, viewport = Viewport.NORTH_WEST)
c2.isometricView()

# Create a map on plane cut instance for the third viewport.
mopc3 = MapOnPlaneCut(scene = s, data_collector = dc2, 
        viewport = Viewport.NORTH_EAST)
mopc3.setPlaneToXY()
mopc3.rotateX(89.9)

c3 = Camera(scene = s, data_collector = dc2, viewport = Viewport.NORTH_EAST)
c3.bottomView()
c3.azimuth(-40)

# Create two map on plance cut instances for the fourth viewport.
mopc4_1 = MapOnPlaneCut(scene = s, data_collector = dc2, 
        viewport = Viewport.SOUTH_EAST)
mopc4_1.setPlaneToXZ()
mopc4_1.rotateZ(-20)
mopc4_1.setOpacity(0.8)

mopc4_2 = MapOnPlaneCut(scene = s, data_collector = dc2, 
        viewport = Viewport.SOUTH_EAST)
mopc4_2.setPlaneToXY()
mopc4_2.rotateY(20)
mopc4_2.setOpacity(0.8)

c4 = Camera(scene = s, data_collector = dc2, viewport = Viewport.SOUTH_EAST)
c4.elevation(-30)

s.render()
