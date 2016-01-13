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

from esys.pyvisi import Scene, Text2D, LocalPosition
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 800, 
        y_size = 600)

# Create one 2D text instance for the first viewport.
t1 = Text2D(scene = s, viewport = Viewport.SOUTH_WEST, text = "VTK ...")
t1.setPosition(LocalPosition(20, 30))
t1.setColor(Color.BLUE)
t1.setFontSize(20)
t1.boldOn()

# Create one 2D text instance for the third viewport.
t2 = Text2D(scene = s, viewport = Viewport.NORTH_EAST, 
        text = "PYTHON ...\n MESH")
t2.setPosition(LocalPosition(20, 30))
t2.setColor(Color.PURPLE)
t2.setFontSize(50)
t2.setFontToArial()
t2.shadowOn()

s.render()

