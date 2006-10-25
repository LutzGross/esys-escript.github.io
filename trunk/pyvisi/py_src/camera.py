"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from geo import Position

class Camera:
	"""
	Class that controls the camera and its settings.
	"""

	def __init__(self, scene, data_collector):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		"""

		self.scene = scene
		self.data_collector = data_collector
		self.vtk_camera = vtk.vtkCamera()		

		self.setCamera()

	def setCamera(self):
		"""
		Set up the camera and associate it with the renderer.
		"""

		self.resetCamera()	
		# Resets the camera clipping range to ensure no objects are cut off.
		self.scene.getRenderer().ResetCameraClippingRange()	
		
		# Default camera focal point and position is the center.
		center = self.data_collector.getReader().GetOutput().GetCenter()
		self.setFocalPoint(Position(center[0], center[1], center[2]))
		#print center[0], center[1], center[2]
		# Camera distance from the rendered object = z + (x*4)
		self.setPosition(
			Position(center[0], center[1], center[2] + (center[0] * 4)))
		# Assigns camera to the scene.
		self.scene.getRenderer().SetActiveCamera(self.vtk_camera)
	
		# Resets the camera clipping range to ensure no objects are cut off.
		self.scene.getRenderer().ResetCameraClippingRange()	
		self.resetCamera()	

	def setClippingRange(self, near_clipping, far_clipping):
		"""
		Set the near and far clipping plane of the camera.
		@type near_clipping: Number
		@param near_clipping: Distance to the near clipping plane
		@type far_clipping: Number
		@param far_clipping: Distance to the far clipping plane
		"""

		self.vtk_camera.SetClippingRange(near_clipping, far_clipping)
		self.resetCamera()

	def setFocalPoint(self, position):
		"""
		Set the focal point of the camera.
		@type position: L{Position <geo.Position>} object
		@param position: Camera focal point 
		"""

		self.vtk_camera.SetFocalPoint(position.getXCoor(), position.getYCoor(),
			position.getZCoor())	
		self.resetCamera()

	def setPosition(self, position):
		"""
		Set the position of the camera.
		@type position: L{Position <geo.Position>} object
		@param position: Camera position 
		"""

		self.vtk_camera.SetPosition(position.getXCoor(), position.getYCoor(),
			position.getZCoor())
		self.resetCamera()

	def setViewUp(self, position):
		"""
		Set the up direction of the camera.
		@type position: L{Position <geo.Position>} object
		@param position: Camera up direction
		"""

		self.vtk_camera.SetViewUp(position.getXCoor(), position.getYCoor(), 
			position.getZCoor())
		self.resetCamera()

	def zoom(self, factor):
		"""
		Zoom in and out of the rendered object.
		@type factor: Number
		@param factor: Amount to zoom in and out
		"""

		self.vtk_camera.Zoom(factor)
		self.resetCamera()

	def azimuth(self, angle):
		"""
		Rotate the camera to the left and right.
		@type angle: Number 
		@param angle: Degree to rotate the camera 
		"""

		self.vtk_camera.Azimuth(angle)
		self.resetCamera()

	def elevation(self, angle):
		"""
		Rotate the camera to the top and bottom.
		@type angle: Number
		@param angle: Degree to rotate the camera
		"""

		self.vtk_camera.Elevation(angle)
		self.resetCamera()


	def roll(self, angle):
		"""
		Roll the camera to the left and right.
		@type angle: Number
		@param angle: Degree to turn the camera
		"""

		self.vtk_camera.Roll(angle)
		self.resetCamera()

	def dolly(self, distance):
		"""
		Move the camera towards and away from the focal point along the normal.
		@type distance: Number
		@param distance: Amount to move towards and away 
		"""
		
		self.vtk_camera(distance)
		self.resetCamera()

	def backView(self):
		"""
		View the back of the rendered object.
		"""
	
		self.azimuth(180)	
		self.resetCamera()

	def topView(self):
		"""
		View the top of the rendered object.
		"""

		self.elevation(90)
		self.resetCamera()

	def bottomView(self):
		"""
		View the bottom of the rendered object.
		"""

		self.elevation(-90)
		self.resetCamera()

	def leftView(self):
		"""
		View the left side of the rendered object.
		"""

		self.azimuth(-90)
		self.resetCamera()

	def rightView(self):
		"""
		View the right side of the rendered object.
		"""

		self.azimuth(90)
		self.resetCamera()

	def resetCamera(self):
		"""
		Repositions the camera so that all actors can be seen. Needs to
		be called whenever the camera's settings are modified.
		"""

		self.scene.getRenderer().ResetCamera()	


	#def isometricView(self):
	#	self.elevation(-50)
		#self.roll(30)
		#self.azimuth(-90)
