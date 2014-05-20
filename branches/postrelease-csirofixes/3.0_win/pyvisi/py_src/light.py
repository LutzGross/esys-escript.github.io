
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
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


from position import GlobalPosition
from constant import Viewport
from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class Light:
	"""
	Class that defines a light. A light controls the lighting for the
	rendered object and works in a similar way to L{Camera <camera.Camera>}.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the light.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to	
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.Light is not running on more than one processor."

		self.__viewport = viewport
		self.__vtk_light = vtk.vtkLight()

		# Keeps track whether light has been modified.
		self.__modified = True
		# Keeps track whether the modification to the light was due to the
		# instantiation. If it is, then __setupLight() method is called.
		self.__initialization = True
		scene._addVisualizationModules(self)

	def __setupLight(self, scene):
		"""
		Set up the light and associate it with the renderer.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		scene._addLight(self.__viewport, self.__vtk_light)

	def setColor(self, color):
		"""
		Set the light color.

		@type color: L{Color <constant.Color>} constant
		@param color: Light color
		"""

		self.__vtk_light.SetColor(color)
		self.__modified = True

	def setFocalPoint(self, position):
		"""
		Set the focal point of the light.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Light focal point
		"""

		self.__vtk_light.SetFocalPoint(position._getGlobalPosition())
		self.__modified = True

	def setPosition(self, position):
		"""
		Set the position of the light.

		@type position: L{GlobalPosition <position.GlobalPosition>} object
		@param position: Light position
		"""

		self.__vtk_light.SetPosition(position._getGlobalPosition())
		self.__modified = True

	# Elevation and azimuth is set to zero so that users do not 
	# have to change both at the same time.
	def setAngle(self, elevation = 0, azimuth = 0):
		"""
		An alternative to set the position and focal point of the light 
		by using the elevation and azimuth.	

		@type elevation: Number
		@param elevation: Degree to rotate the light to the top and bottom
		@type azimuth: Number
		@param azimuth: Degree to rotate the light to the left and right
		"""

		# NOTE: The elevation angle of light does not appear to suffer the same
		# constraint as the elevation angle of camera where the elevation
		# angle is constraint between 90/-90.
		self.__vtk_light.SetDirectionAngle(elevation, azimuth)
		self.__modified = True

	def setIntensity(self, intensity):
		"""
		Set the intensity (brightness) of the light.

		@type intensity: Number
		@param intensity: Intensity (brightness) of the light
		"""

		self.__vtk_light.SetIntensity(intensity)
		self.__modified = True

	def _isModified(self):
		"""
		Return whether the Light has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		if (self.__modified == True):
			return True
		else:
			return False

	def _render(self, scene):
		"""
		Render the light.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if(self._isModified() == True):
			# Will only be true once only when the light is instantiated.
			if(self.__initialization == True): 
				self.__setupLight(scene)
				self.__initialization == False

			self.__modified = False
			

