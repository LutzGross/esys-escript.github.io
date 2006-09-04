"""
Class and functions associated with cameras and views
"""

import vtk

class Camera:
	def __init__(self, open_scene):
		self.open_scene = open_scene
		self.vtk_camera = None 

		self.setCamera()

	def setCamera(self):
		self.vtk_camera = vtk.vtkCamera()		
		self.open_scene.getRenderer().SetActiveCamera(self.vtk_camera)

	def setClippingRange(self, near_clipping, far_clipping):
		self.vtk_camera.SetClippingRange(near_clipping, far_clipping)

	def setFocalPoint(self, position):
		self.vtk_camera.SetFocalPoint(position.getXCoor(), position.getYCoor(),
			position.getZCoor())	

	def setPosition(self, position):
		self.vtk_camera.SetPosition(position.getXCoor(), position.getYCoor(),
			position.getZCoor())

	def setViewUp(self, position):
		self.vtk_camera.SetViewUp(position.getXCoor(), position.getYCoor(), 
			position.getZCoor())






class FrontView(Camera):
    pass

class BackView(Camera):
    pass

class TopView(Camera):
    pass

class BottomView(Camera):
    pass

class LeftView(Camera):
    pass

class RightView(Camera):
    pass

class IsometricView(Camera):
    pass
