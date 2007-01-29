"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from actor import Actor3D
from constant import Viewport
from texture import Texture
from plane import PlaneSource

# NOTE: DataSetMapper, Actor3D, Texture and PlaneSource were inherited to 
# allow access to their public methods from the driver.
class Image(DataSetMapper, Actor3D, Texture, PlaneSource):
	"""
	Class that displayes an image with interaction capability.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, image_reader, viewport = Viewport.SOUTH_WEST):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type image_reader: L{ImageReader <datacollector.ImageReader>}
				object
		@param image_reader: Deal with source of image for visualisation
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which the object is to be rendered on 
		"""

		# ----- Image -----
		print "image"
		Texture.__init__(self, image_reader._getOutput())
		PlaneSource.__init__(self)
		DataSetMapper.__init__(self, PlaneSource._getOutput(self))	

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		Actor3D._setTexture(self, Texture._getTexture(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))
		


