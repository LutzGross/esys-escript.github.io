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
from actor import Actor2D
from constant import Viewport, Color

# NOTE: Actor2D was inherited to allow access to its public methods from 
# the driver.
class Text2D(Actor2D):
	"""
	Class that defines a 2D text actor. A two-dimensional text is used to
	annotate the rendered object (i.e. inserting titles, authors and labels).
	"""

	def __init__(self, scene, text, viewport = Viewport.SOUTH_WEST):
		"""
		Initialise the 2D text actor.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type text: String
		@param text: 2D text to be displayed
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		"""

		self.__text = text
		self.__viewport = viewport
		self._vtk_actor2D = vtk.vtkTextActor()

		self.__setupText2D(scene)
	
	def __setupText2D(self, scene):
		"""
		Setup the 2D text.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		self.__setInput()
		# Set the color of the 2D text to black by default.
		self.setColor(Color.BLACK)
		# Add the 2D text to the appropriate renderer.
		scene._addActor2D(self.__viewport, self._vtk_actor2D)

	def __setInput(self):
		"""
		Set the input for the 2D text.
		"""

		self._vtk_actor2D.SetInput(self.__text)

	def setFontSize(self, size):
		"""
		Set the 2D text size.

		@type size: Number
		@param size: Size of the 2D text
		"""

		self._vtk_actor2D.GetTextProperty().SetFontSize(size)

	def setFontToTimes(self):
		"""
		Set the 2D text font type to Times New Roman.
		"""

		self._vtk_actor2D.GetTextProperty().SetFontFamilyToTimes()

	def setFontToArial(self):
		"""
		Set the 2D text font type to Arial.
		"""

		self._vtk_actor2D.GetTextProperty().SetFontFamilyToArial()

	def setFontToCourier(self):
		"""
		Set the 2D text front type to Courier.
		"""

		self._vtk_actor2D.GetTextProperty().SetFontFamilyToCourier()

	def boldOn(self):
		"""
		Bold the 2D text.
		"""

		self._vtk_actor2D.GetTextProperty().BoldOn()

	def shadowOn(self):
		"""
		Apply shadow onto the 2D text to ease visibility.
		"""

		self._vtk_actor2D.GetTextProperty().ShadowOn()

	def setColor(self, color):
		"""
		Set the color of the 2D text.

		@type color: L{Color <constant.Color>} constant
		@param color: 2D text color
		"""

		self._vtk_actor2D.GetTextProperty().SetColor(color)

