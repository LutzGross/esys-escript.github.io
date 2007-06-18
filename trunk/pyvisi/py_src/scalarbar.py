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
from position import LocalPosition
from constant import Color

class ScalarBar:
	"""
	Class that defines a scalar bar.
	"""
	
	def __init__(self):
		"""
		Initialise the scalar bar.
		"""

		self.__vtk_scalar_bar = vtk.vtkScalarBarActor()
		self.__setupScalarBar()

	def __setupScalarBar(self):
		"""
		Setup the scalar bar.
		"""

		self.__vtk_scalar_bar.GetPositionCoordinate().SetCoordinateSystemToViewport()
		self.setNumberOfLabels(8)
		# Set the default color of the scalar bar's title and label to black.
		self.setTitleColor(Color.BLACK)
		self.setLabelColor(Color.BLACK)
		self.titleBoldOff()
		self.labelBoldOff()

	def _setScalarBarLookupTable(self, lookup_table):	
		"""
		Set the scalar bar's lookup table.

		@type lookup_table: vtkScalarsToColors
		@param lookup_table: Converts scalar data to colors
		"""

		self.__vtk_scalar_bar.SetLookupTable(lookup_table)

	def setNumberOfLabels(self, labels):
		"""
		Set the number of lables on the scalar bar.

		@type labels: Number
		@param labels: Number of labels on the scalar bar
		"""

		self.__vtk_scalar_bar.SetNumberOfLabels(labels)

	def setTitle(self, title):	
		"""
		Set the title on the scalar bar.

		@type title: String
		@param title: Title on the scalar bar
		"""

		self.__vtk_scalar_bar.SetTitle(title)

	def setPosition(self, position):
		"""
		Set the local position of the scalar bar.

		@type position: L{LocalPosition <position.LocalPosition>} object
		@param position: Position of the scalar bar
		"""

		self.__vtk_scalar_bar.SetPosition(position._getXCoor(), \
				position._getYCoor())

	def setOrientationToHorizontal(self):
		"""
		Set the orientation of the scalar bar to horizontal.
		"""

		self.__vtk_scalar_bar.SetOrientationToHorizontal()
		self.setWidth(0.8)
		self.setHeight(0.10)
		self.setPosition(LocalPosition(130,5))

	def setOrientationToVertical(self):
		"""
		Set the orientation of the scalar bar to vertical.
		"""

		self.__vtk_scalar_bar.SetOrientationToVertical()
		self.setWidth(0.14)
		self.setHeight(0.85)
		self.setPosition(LocalPosition(5,150))

	def setHeight(self, height):
		"""
		Set the height of the scalar bar.
	
		@type height: Number
		@param height: Height of the scalar bar
		"""

		self.__vtk_scalar_bar.SetHeight(height)

	def setWidth(self, width):
		"""
		Set the width of the scalar bar.

		@type width: Number
		@param width: Width of the scalar bar.
		"""

		self.__vtk_scalar_bar.SetWidth(width)

	def setLabelColor(self, color):
		"""
		Set the color of the scalar bar's label.

		@type color: L{Color <constant.Color>} constant
		@param color: Color of the scalar bar label
		"""

		self.__getLabelProperty().SetColor(color)

	def labelBoldOn(self):
		"""
		Enable the bold on the scalar bar's label.
		"""

		self.__getLabelProperty().BoldOn()

	def labelBoldOff(self):
		"""
		Disable the bold on the scalar bar's label.
		"""

		self.__getLabelProperty().BoldOff()

	def setTitleColor(self, color):
		"""
		Set the color of the scalar bar's title.

		@type color: L{Color <constant.Color>} constant
		@param color: Color of the scalar bar title
		"""

		self.__getTitleProperty().SetColor(color)

	def titleBoldOn(self):
		"""
		Enable the bold on the scalar bar's title.
		"""

		self.__getTitleProperty().BoldOn()

	def titleBoldOff(self):
		"""
		Disable the bold on the scalar bar's title.
		"""

		self.__getTitleProperty().BoldOff()

	def __getLabelProperty(self):
		"""
		Return the scalar bar's label property.

		@rtype: vtkTextPorperty
		@return: Scalar bar's label property 
		"""

		return self.__vtk_scalar_bar.GetLabelTextProperty()

	def __getTitleProperty(self):
		"""
		Return the scalar bar's title property.

		@rtype: vtkTextPorperty
		@return: Scalar bar's title property 
		"""

		return self.__vtk_scalar_bar.GetTitleTextProperty()

	def _getScalarBar(self):
		"""
		Return the scalar bar.

		@rtype: vtkScalarBarActor
		@return: Scalar bar actor
		"""

		return self.__vtk_scalar_bar
	
