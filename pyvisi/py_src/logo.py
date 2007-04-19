"""
@author: John NGUI
"""

import vtk
from mapper import ImageMapper
from imagereslice import ImageReslice
from actor import Actor2D
from constant import Viewport

# NOTE: ImageMapper, ImageReslice and Actor2D were inherited to  allow access 
# to their public methods from the driver.
class Logo(ImageMapper, ImageReslice, Actor2D):
	"""
	Class that displays an image which can be scaled (upwards and downwards). 
	The image can also be translated and rotated along the X, Y and Z axes.

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
		@param scene: Scene in which objects are to be rendered on
		@type image_reader: L{ImageReader <imagereader.ImageReader>}
				object
		@param image_reader: Deal with source of image for vizualisation
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which objects are to be rendered on 
		"""

		# ----- Image -----

		ImageReslice.__init__(self, image_reader._getOutput())
		ImageMapper.__init__(self, ImageReslice._getOutput(self))

		Actor2D.__init__(self, ImageMapper._getImageMapper(self))
		scene._addActor2D(viewport, Actor2D._getActor2D(self))

