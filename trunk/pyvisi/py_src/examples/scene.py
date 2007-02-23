from esys.pyvisi import Scene
from esys.pyvisi.constant import *

# Create a scene instance.
s = Scene(renderer = Renderer.ONLINE, num_viewport = 4, x_size = 800, 
        y_size = 600)
s.setBackground(color = Color.YELLOW)
s.setTitleBar(text = "Prototype...")

# Render the scene.
s.render()

