
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
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut, ColorMode 
from streamlinemodule import  StreamLineModule
from tube import Tube
from point import PointSource
from outline import Outline
from average import CellDataToPointData
from position import GlobalPosition
from esys.escript import getMPISizeWorld

# NOTE: DataSetMapper, Actor3D, PointSource, StreamLineModule and Tube  were 
# inherited to allow access to their public methods from the driver.
class StreamLine(DataSetMapper, Actor3D, PointSource, StreamLineModule, Tube):
	"""
	Class that shows the direction of particles of a vector field using
	streamlines.The streamlines can either be colored or grey-scaled,
	depending on the lookup table used. If the streamlines are colored,
	there are two possible coloring modes: (1) using vector data or (2)
	using scalar data.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
		color_mode = ColorMode.VECTOR, lut = Lut.COLOR, cell_to_point = False, 
		outline = True): 
		"""
		Initialise the StreamLine.

		@attention: The source can either be point or cell data. If the 
		source is cell data, a conversion to point data may or may not be 
		required, in order for the object to be rendered correctly. 
		If a conversion is needed, the 'cell_to_point' flag must be set to 
		'True', otherwise 'False' (which is the default). On occasions, an
		inaccurate object may be rendered from cell data even after conversion.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__color_mode = color_mode
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline

		# Keeps track whether Streamline has been modified.
		self.__modified = True 
		PointSource.__init__(self)
		StreamLineModule.__init__(self)
		Tube.__init__(self)
		DataSetMapper.__init__(self)
		Actor3D.__init__(self)
		scene._addVisualizationModules(self)

		# ----- Outline -----

		# NOTE: Changes cannot be made to the Outline's properties from the 
		# driver.
		if(self.__outline == True):
			outline = Outline(self.__data_collector._getDataCollectorOutput())
			mapper = DataSetMapper()
			mapper._setupDataSetMapper(outline._getOutlineOutput()) 

			actor3D = Actor3D()
			actor3D._setupActor3D(mapper._getDataSetMapper())
			# Default outline color is black.
			actor3D.setColor(Color.BLACK)

			# Default line width is 1.
			actor3D._setLineWidth(1)
			scene._addActor3D(self.__viewport, actor3D._getActor3D())

		# ----- Streamline -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(self.__lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(self.__lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		if(self.__cell_to_point == True): # Converts cell data to point data.
			c2p = CellDataToPointData(
					self.__data_collector._getDataCollectorOutput())
			self._setupPointSource(c2p._getCellToPointOutput())
			self._setupStreamLineModule(c2p._getCellToPointOutput(), 
					self._getPointSourceOutput())
		elif(self.__cell_to_point == False): # No conversion happens.
			self._setupPointSource(
					self.__data_collector._getDataCollectorOutput())
			self._setupStreamLineModule(
					self.__data_collector._getDataCollectorOutput(), 
					self._getPointSourceOutput())

		self._setupTube(self._getStreamLineModuleOutput())
		self._setupDataSetMapper(self._getTubeOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the StreamLine or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the streamline.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self._isPointSourceCenterSet() != True):
				center = self.__data_collector._getCenter()
				center = GlobalPosition(center[0], center[1], center[2])
				self.setPointSourceCenter(center)

			self._setPointSourceCenter()	

			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			if(self.__data_collector._isVectorSet() == True):
				self.__data_collector._setActiveVector()

			# Color streamline by vector.
			if(self.__color_mode == ColorMode.VECTOR): 				
				self._setScalarVisibilityOn()
				self._setSpeedScalarsOn()

				# self._isScalarRangeSet checks whether the scalar range has 
				# beenspecified by the user. If it has, then the scalar range
				# read from the source will be ignored.
				if(not(self._isScalarRangeSet())): 
					self._setScalarRange(\
							self.__data_collector._getVectorRange())
			# Color streamline by scalar.
			elif(self.__color_mode == ColorMode.SCALAR): 				
				self._setScalarVisibilityOn()
				self._setSpeedScalarsOff()

				# self._isScalarRangeSet checks whether the scalar range has 
				# beenspecified by the user. If it has, then the scalar range
				# read from the source will be ignored.
				if(not(self._isScalarRangeSet())): 
					self._setScalarRange(\
							self.__data_collector._getScalarRange())

			self.__modified = False


