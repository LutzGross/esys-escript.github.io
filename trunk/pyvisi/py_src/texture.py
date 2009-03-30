
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
from esys.escript import getMPISizeWorld

class Texture:
	"""
	Class that defines the texture for the rendered object.
	"""

	def __init__(self):
		"""
		Initialise the texture.
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Texture is not running on more than one processor."
		self.__vtk_texture = vtk.vtkTexture()

	def _setupTexture(self, image):
		"""
		Setup the texture.

		@type image: vtkImageData
		@param image: Image from which data is to be read
		"""

		self.__image = image
		self.__setInput()	

	def __setInput(self):
		"""
		Set the input for the texture.
		"""

		self.__vtk_texture.SetInput(self.__image)

	def _getTexture(self):
		"""
		Return the texture.

		@rtype: vtkTexture
		@return: Texture of the rendered object
		"""

		return self.__vtk_texture

