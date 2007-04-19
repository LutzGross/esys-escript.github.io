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

		# ----- Image -----

		Texture.__init__(self, image_reader._getOutput())
		PlaneSource.__init__(self)

		Transform.__init__(self)
		TransformFilter.__init__(self, PlaneSource._getOutput(self),
				Transform._getTransform(self))

		DataSetMapper.__init__(self, TransformFilter._getOutput(self))

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		Actor3D._setTexture(self, Texture._getTexture(self))

		scene._addActor3D(viewport, Actor3D._getActor3D(self))

