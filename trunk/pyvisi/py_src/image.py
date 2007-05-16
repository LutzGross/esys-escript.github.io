"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from actor import Actor3D
from constant import Viewport
from texture import Texture
from plane import PlaneSource
from transform import Transform, TransformFilter

# NOTE: DataSetMapper, Actor3D, Texture, PlaneSource, Transform and
# TransformFilter were inherited to  allow access to their public methods 
# from the driver.
class Image(DataSetMapper, Actor3D, Texture, PlaneSource, Transform,
		TransformFilter):
	"""
	Class that displays an image which can be scaled (upwards and downwards)
	and has interaction capability. The image can also be translated and 
	rotated along the X, Y and Z axes.

	@bug: Translating an image works differently (opposite) compared to 
	translating a plane. For example, a positive translation along the 
	z-axis moves a plane up. However, if the identical translation is applied on
	an image, the image moves down.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, image_reader, viewport = Viewport.SOUTH_WEST):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which the image is to be displayed
		@type image_reader: L{ImageReader <imagereader.ImageReader>}
				object
		@param image_reader: Deal with source of data for vizualisation
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which the image is to be displayed
		"""

		self.__scene = scene
		self.__image_reader = image_reader
		self.__viewport = viewport

		# Keeps track whether Image has been modified.
		self.__modified = True 
		Texture.__init__(self)
		PlaneSource.__init__(self)
		Transform.__init__(self)
		TransformFilter.__init__(self)
		DataSetMapper.__init__(self)
		Actor3D.__init__(self)
		scene._addVisualizationModules(self)

		# ----- Image -----

		self._setupTexture(image_reader._getImageReaderOutput())
		self._setupTransformFilter(self._getPlaneSourceOutput(),
				self._getTransform())

		self._setupDataSetMapper(self._getTransformFilterOutput())
		self._setupActor3D(self._getDataSetMapper())

		self._setTexture(self._getTexture())
		self.__scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):
		"""
		Return whether the Image has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		if (self.__modified == True):
			return True
		else:
			return False

	def _render(self):
		"""
		Render the image.
		"""

		if(self._isModified() == True):
			self.__isModified = False
			

