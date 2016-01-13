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

from esys.pyvisi import Scene
from esys.pyvisi.constant import *

# Create a scene instance.
s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 800, 
        y_size = 600)
s.setBackground(Color.YELLOW)
s.setTitleBar("Prototype...")

# Render the scene.
s.render()

