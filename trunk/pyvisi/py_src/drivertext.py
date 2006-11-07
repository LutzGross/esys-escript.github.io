from scene import Scene
from text import Text
from constants import *
from style import Style

s = Scene(renderer = "vtk_online", x_size = 800, y_size = 600)

# Set the style of the text.
st = Style()
st.setColor(BLUE)
st.setFontFamily("Arial")
st.boldOn()
st.shadowOn()
st.italicOn()

# Create a 2D text.
t = Text(scene = s)
t.setText("VTK Rendering...")
t.setPosition(120, 300)
t.setStyle(st)
s.render()

