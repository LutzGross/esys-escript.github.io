
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


from mapper import DataSetMapper
from actor import Actor3D
from constant import Viewport
from cube import CubeSource
from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

# NOTE: CubeSource, DataSetMapper and Actor3D were inherited to allow 
# access to their public methods from the driver.
class Rectangle(CubeSource, DataSetMapper, Actor3D):
	"""
	Class that generates a rectangle box.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the Rectangle.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which objects are to be rendered on 
		"""

		self.__viewport = viewport
		
		CubeSource.__init__(self)	
		DataSetMapper.__init__(self)
		Actor3D.__init__(self)

		# ----- Rectangle -----

		self._setupDataSetMapper(self._getCubeSourceOutput())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())
	


