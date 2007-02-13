from esys.pyvisi import Scene
from esys.pyvisi.constant import *

# Create a scene instance.
s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 800, 
        y_size = 600)
s.setBackground(Color.YELLOW)
s.setTitleBar("Prototype...")

# Render the scene.
s.render()

