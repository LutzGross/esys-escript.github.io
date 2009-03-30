
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
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


import vtk
from position import GlobalPosition
from constant import Viewport
from esys.escript import getMPISizeWorld

class Camera:
	"""
	Class that defines a camera. A camera controls the display angle of
	the rendered object and one is usually created for a
	L{Scene <scene.Scene>}. However, if a L{Scene <scene.Scene>} has four
	viewports, then a separate camera may be created for each viewport.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the camera.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Camera is not running on more than one processor."
		self.__viewport = viewport
		self.__vtk_camera = vtk.vtkCamera()

		# Keeps track whether camera has been modified.
		self.__modified = True 
		# Keeps track whether the modification to the camera was due to the 
		# instantiation. If it is, then __setupCamera() method is called.
		self.__initialization = True
		scene._addVisualizationModules(self)

	def __setupCamera(self, scene):
		"""
		Setup the camera.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		# Assign the camera to the appropriate renderer
		scene._setActiveCamera(self.__viewport, self.__vtk_camera)
		self.__resetCamera(scene)	

	def setFocalPoint(self, position):
		"""
		Set the focal point of the camera.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera focal point
		"""

		self.__vtk_camera.SetFocalPoint(position._getGlobalPosition())		
		self.__modified = True 

	def setPosition(self, position):
		"""
		Set the position of the camera.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera position
		"""

		self.__vtk_camera.SetPosition(position._getGlobalPosition())
		self.__modified = True 

	def setClippingRange(self, near_clipping, far_clipping):
		"""
		Set the near and far clipping plane of the camera.

		@type near_clipping: Number
		@param near_clipping: Distance to the near clipping plane
		@type far_clipping: Number
		@param far_clipping: Distance to the far clipping plane
		"""
	
		self.__vtk_camera.SetClippingRange(near_clipping, far_clipping)

	def setViewUp(self, position):
		"""
		Set the view up direction of the camera.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Camera view up direction
		"""

		self.__vtk_camera.SetViewUp(position._getGlobalPosition())

	def azimuth(self, angle):
		"""
		Rotate the camera to the left and right.

		@type angle: Number
		@param angle: Degree to rotate the camera
		"""

		self.__vtk_camera.Azimuth(angle)
		self.__modified = True 

	def elevation(self, angle):
		"""
		Rotate the camera to the top and bottom.
		
		@type angle: Number
		@param angle: Degree to rotate the camera (only between -90 and 90)
		"""

		self.__vtk_camera.Elevation(angle)
		# Recompute the view up vector. If not used the elevation angle is  
		# unable to exceed 87/-87 degrees. A warning resetting the
		# view up will also be thrown and the rendered object may be incorrect.
		# With the view up recomputed, the elevation angle can reach between
		# 90/-90 degrees. Exceeding that, the rendered object may be incorrect.
		self.__vtk_camera.OrthogonalizeViewUp()
		self.__modified = True 

	def roll(self, angle):
		"""
		Roll the camera to the left and right.

		@type angle: Number
		@param angle: Degree to roll the camera
		"""

		self.__vtk_camera.Roll(-angle)
		self.__modified = True 

	def backView(self):
		"""
		Rotate the camera to view the back of the rendered object.
		"""

		self.azimuth(180)

	def topView(self):
		"""
		Rotate the camera to view the top of the rendered object.
		"""
		
		self.elevation(90)

	def bottomView(self):
		"""
		Rotate the camera to view the bottom of the rendered object.
		"""

		self.elevation(-90)

	def leftView(self):
		"""
		Rotate the camera to view the left side of the rendered object.
		"""

		self.azimuth(-90)

	def rightView(self):
		"""
		Rotate the camera to view the right side of the rendered object.
		"""

		self.azimuth(90)

	def isometricView(self):
		"""
		Rotate the camera to view the isometric angle of the rendered object.
		"""

		self.roll(-45)
		self.elevation(-45)
		
	def dolly(self, distance):
		"""
		Move the camera towards (greater than 1) the rendered object. However,
		the camera is unable to be moved away from the rendered object.

		@type distance: Number
		@param distance: Amount to move towards the rendered object
		"""

		self.__vtk_camera.Dolly(distance)
		self.__modified = True 

	def parallelProjectionOn(self):
		"""
		Enable camera parallel projection.
		"""

		self.__vtk_camera.ParallelProjectionOn()
		self.__modified = True 

	def parallelProjectionOff(self):
		"""
		Disable camera parallel projection.
		"""

		self.__vtk_camera.ParallelProjectionOff()
		self.__modified = True

	def __resetCameraClippingRange(self, scene):
		"""
		Reset the camera clipping range based on the bounds of the visible 
		actors. This ensures the rendered object is not cut-off.
		Needs to be called whenever the camera's settings are modified.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		scene._getRenderer()[self.__viewport].ResetCameraClippingRange() 

	def __resetCamera(self, scene):
		"""
		Repositions the camera to view the center point of the actors.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		scene._getRenderer()[self.__viewport].ResetCamera() 

	def _isModified(self):
		"""
		Return whether the Camera has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		if (self.__modified == True):
			return True
		else:
			return False

	def _render(self, scene):
		"""
		Render the camera.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if(self._isModified() == True):
			# Will only be true once only when the camera is instantiated.
			if(self.__initialization == True): 
				self.__setupCamera(scene)
				self.__initialization == False

			self.__resetCameraClippingRange(scene)
			self.__modified = False
			

