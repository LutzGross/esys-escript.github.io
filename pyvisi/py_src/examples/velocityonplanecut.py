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

from esys.pyvisi import Scene, DataCollector, VelocityOnPlaneCut
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = "/home/jongui/data/laurent/slab/source/slab.xml")

# NOTE: There is a difference between performing rotation then followed by
# translation and performing translation then followed by rotation.

# Create a velocity instance.
vopc1 = VelocityOnPlaneCut(scene = s, data_collector = dc1, 
        color_mode = ColorMode.VECTOR)
vopc1.setScaleFactor(2000000)
vopc1.setPlaneToXY(0.2)
vopc1.setDimension(4,4,4)

s.render()

