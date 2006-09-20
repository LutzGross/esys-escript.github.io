"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk

class Camera:
	"""
	Class that controls the camera and its settings.
	"""

	def __init__(self, scene):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		"""

		self.scene = scene
		self.vtk_camera = vtk.vtkCamera()		

		self.setCamera()

	def setCamera(self):
		"""
		Set up the camera and associate it with the renderer.
		"""

		self.scene.getRenderer().SetActiveCamera(self.vtk_camera)

	def setClippingRange(self, near_clipping, far_clipping):
		"""
		Set the near and far clipping plane of the camera.

		@type near_clipping: Number
		@param near_clipping: Distance to the near clipping range
		@type far_clipping: Number
		@param far_clipping: Distance to the far clipping plane
		"""

		self.vtk_camera.SetClippingRange(near_clipping, far_clipping)

	def setFocalPoint(self, position):
		"""
		Set the focal point of the camera.
	
		@type position: L{Position <geo.Position>} object
		@param position: Camera focal point position 
		"""

		self.vtk_camera.SetFocalPoint(position.getXCoor(), position.getYCoor(),
			position.getZCoor())	

	def setPosition(self, position):
		"""
		Set the position of the camera.

		@type position: L{Position <geo.Position>} object
		@param position: Camera position 
		"""

		self.vtk_camera.SetPosition(position.getXCoor(), position.getYCoor(),
			position.getZCoor())

	def setViewUp(self, position):
		"""
		Set the up direction of the camera.

		@type position: L{Position <geo.Position>} object
		@param position: Camera view up position
		"""

		self.vtk_camera.SetViewUp(position.getXCoor(), position.getYCoor(), 
			position.getZCoor())

	def setZoom(self, zoom):
		self.vtk_camera.Zoom(zoom)



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
