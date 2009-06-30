
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Lutz Gross, l.gross@uq.edu.au
Author: John Ngui, john.ngui@uq.edu.au
Author: Ken Steube
"""

# Import the necesary modules.
from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.finley import Rectangle
from esys.pyvisi import Scene, DataCollector, Map, Camera
from esys.pyvisi.constant import *
import os

PYVISI_EXAMPLE_IMAGES_PATH = "images_out"
if not os.path.isdir(PYVISI_EXAMPLE_IMAGES_PATH) and getMPIRankWorld()==0: os.mkdir(PYVISI_EXAMPLE_IMAGES_PATH)

X_SIZE = 400
Y_SIZE = 300
JPG_RENDERER = Renderer.OFFLINE_JPG # change to Renderer.ONLINE_JPG to interact with visualiztion window

# generate domain:
mydomain = Rectangle(l0=1.,l1=1.,n0=40, n1=20)
# define characteristic function of Gamma
x = mydomain.getX()
gammaD = whereZero(x[0])+whereZero(x[1])
# define PDE and get its solution u
mypde = Poisson(domain=mydomain)
mypde.setValue(f=1,q=gammaD)
u = mypde.getSolution()

# Create a Scene.
s = Scene(renderer = JPG_RENDERER, x_size = X_SIZE, y_size = Y_SIZE)

# Create a DataCollector reading directly from an escript object.
dc = DataCollector(source = Source.ESCRIPT)
dc.setData(sol = u)

# Create a Map.
Map(scene = s, data_collector = dc, viewport = Viewport.SOUTH_WEST, 
	  lut = Lut.COLOR, cell_to_point = False, outline = True)

# Create a Camera.
c = Camera(scene = s, viewport = Viewport.SOUTH_WEST)

# Render the object.
s.render(image_name = "offscreen.jpg")
