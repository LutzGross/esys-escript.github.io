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
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut, WarpMode, ColorMode
from warp import Warp
from outline import Outline
from transform import Transform
from plane import Plane
from cutter import Cutter
from average import CellDataToPointData

# NOTE: DataSetMapper, Actor3D, Warp, Transform, Plane and Cutter  were 
# inherited to allow access to their public methods from the driver.
class Carpet(DataSetMapper, Actor3D, Warp, Transform, Plane, Cutter):
	"""
	Class that shows a scalar field on a plane deformated along the normal.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no warp_mode is specified, the data will be deformated using scalar
	# data. If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			warp_mode = WarpMode.SCALAR, lut = Lut.COLOR, 
			cell_to_point = False, outline = True):
		"""
		Initialise the Carpet.

		@attention: The source can either be point or cell data. If the 
		source is cell data, a conversion to point data may or may not be 
		required, in order for the object to be rendered correctly. 
		If a conversion is needed, the 'cell_to_point' flag must be set to 
		'True', otherwise 'False' (which is the default).

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		@param warp_mode: L{WarpMode <constant.WarpMode>} constant
		@type warp_mode: Mode in which to deform the scalar field 
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__warp_mode = warp_mode
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		# Keeps track whether Carpet has been modified.
		self.__modified = True 
		Transform.__init__(self)	
		Plane.__init__(self)
		Cutter.__init__(self)
		Warp.__init__(self, self.__warp_mode)
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

		# ----- Carpet -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(self.__lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(self.__lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		self._setupPlane(self._getTransform())

		if(self.__cell_to_point == True): # Converts cell data to point data.
			c2p = CellDataToPointData(\
					self.__data_collector._getDataCollectorOutput())
			self._setupCutter(c2p._getCellToPointOutput(), self._getPlane()) 	
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupCutter(self.__data_collector._getDataCollectorOutput(), 
					self._getPlane()) 	

		self._setupWarp(self._getCutterOutput())
		self._setupDataSetMapper(self._getWarpOutput(),
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the Carpet or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the carpet.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()

			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
				self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


