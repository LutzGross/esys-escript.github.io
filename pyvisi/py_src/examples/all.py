from esys.pyvisi import Scene, DataCollector, MapOnPlaneCut, Camera 
from esys.pyvisi import VelocityOnPlaneCut, StreamLine, EllipsoidOnPlaneCut
from esys.pyvisi import ContourOnPlaneClip, Text2D, LocalPosition
from esys.pyvisi.constant import *

s = Scene(renderer = Renderer.ONLINE, num_viewport = 1, x_size = 1152, 
        y_size = 864)

dc1 = DataCollector(source = Source.XML)
dc1.setFileName(file_name = 
        "/home/jongui/trunk/pyvisi/test/python/data_data/interior_3D.xml")

mopc1 = MapOnPlaneCut(scene = s, data_collector = dc1, 
        viewport = Viewport.SOUTH_WEST)
mopc1.setPlaneToXY()

vopc1 = VelocityOnPlaneCut(scene = s, data_collector = dc1,
        arrow = Arrow.THREE_D, color_mode = ColorMode.SCALAR)
vopc1.setScaleFactor(0.2)
vopc1.setPlaneToYZ(2.999)
vopc1.setDimension(1.5,1.5,1.5)

sl1 = StreamLine(scene = s, data_collector = dc1,
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True,
		color_mode = ColorMode.SCALAR)
sl1.setTubeRadius(0.02)

eopc1 = EllipsoidOnPlaneCut(scene = s, data_collector = dc1, tensor = None,
        viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
eopc1.setScaleFactor(0.1)
eopc1.setPlaneToXZ()
eopc1.rotateX(-45)
eopc1.translate(0,1.6,0)
eopc1.setDimension(1,1,1)

ctropc1 = ContourOnPlaneClip(scene = s, data_collector = dc1, scalar = None,
        viewport  = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True)
ctropc1.setPlaneToXY()
ctropc1.rotateY(10)
ctropc1.generateContours(3)

t1 = Text2D(scene = s, viewport = Viewport.SOUTH_WEST, text = "Pyvisi")
t1.setPosition(LocalPosition(530, 30))
t1.setColor(Color.BLACK)
t1.setFontSize(30)
t1.boldOn()

c1 = Camera(scene = s, data_collector = dc1, viewport = Viewport.SOUTH_WEST)
c1.isometricView()

s.render()
