
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
from mapper import ImageMapper
from imagereslice import ImageReslice
from actor import Actor2D
from constant import Viewport
from esys.escript import getMPISizeWorld

# NOTE: ImageMapper, ImageReslice and Actor2D were inherited to allow access 
# to their public methods from the driver.
class Logo(ImageMapper, ImageReslice, Actor2D):
	"""
	Class that displays a static image, in particular a logo 
	(i.e. company symbol) and has NO interaction capability. The position and
	size of the logo can be specified.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, image_reader, viewport = Viewport.SOUTH_WEST):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which the logo is to be displayed
		@type image_reader: L{ImageReader <imagereader.ImageReader>}
		object
		@param image_reader: Deal with source of data for vizualisation
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which the logo is to be displayed
		"""

		self.__image_reader = image_reader
		self.__viewport = viewport

		self.__modified = True # Keeps track whether Logo has been modified.
		ImageReslice.__init__(self)
		ImageMapper.__init__(self)
		Actor2D.__init__(self)
		scene._addVisualizationModules(self)

		# ----- Logo -----

		self._setupImageReslice(self.__image_reader._getImageReaderOutput())
		self._setupImageMapper(self._getImageResliceOutput())

		self._setupActor2D(self._getImageMapper())
		scene._addActor2D(self.__viewport, self._getActor2D())

	def _isModified(self):	
		"""
		Return whether the Logo or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the logo.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which the logo is to be displayed
		"""

		if (self._isModified() == True):
			self.__modified = False



