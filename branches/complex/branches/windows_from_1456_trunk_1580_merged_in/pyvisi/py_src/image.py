"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


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

	@attention: Translating an image works differently (opposite) compared to 
	translating a plane. For example, a positive translation along the 
	z-axis moves a plane up. However, if the identical translation is applied 
	on an image, the image moves down.
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
		scene._addActor3D(self.__viewport, self._getActor3D())

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

	def _render(self, scene):
		"""
		Render the image.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which the image is to be displayed
		"""

		if(self._isModified() == True):
			self.__isModified = False
			

