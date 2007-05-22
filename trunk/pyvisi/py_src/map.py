"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from actor import Actor3D
from lookuptable import LookupTable
from outline import Outline
from constant import Viewport, Color, Lut, ColorMode
from average import CellDataToPointData

# NOTE: DataSetMapper and Actor3D were inherited to allow access to their 
# public methods from the driver.
class Map(DataSetMapper, Actor3D):
	"""
	Class that shows a scalar field on a domain surface. The domain surface
	can either be colored or grey-scaled, depending on the lookup table used.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):
		"""
		Initialise the Map.

		@attention: The source can either be point or cell data. If the 
		source is cell data, a conversion to point data may or may not be 
		required, in order for the object to be rendered correctly. 
		If a conversion is needed, the 'cell_to_point' flag must be set to 
		'True', otherwise 'False' (which is the default).

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for vizualisation
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which objects are to be rendered on 
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		self.__modified = True # Keeps track whether Map has been modified.
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

		# ----- Map -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(self.__lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable() 
			lookup_table._setTableValue()
		elif(self.__lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable() 
			lookup_table._setLookupTableToGreyScale()

		if(self.__cell_to_point == True): # Convert cell data to point data.
			c2p = CellDataToPointData(
					self.__data_collector._getDataCollectorOutput())
			self._setupDataSetMapper(c2p._getCellToPointOutput(), 
					lookup_table._getLookupTable())	
		elif(self.__cell_to_point == False): # No conversion happens.
			self._setupDataSetMapper(
					self.__data_collector._getDataCollectorOutput(), 
					lookup_table._getLookupTable())	

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())
	
	def _isModified(self):	
		"""
		Return whether the Map or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the surface map.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


###############################################################################


from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, Transform, Plane and Cutter were inherited 
# to allow access to their public methods from the driver.
class MapOnPlaneCut(DataSetMapper, Actor3D, Transform, Plane, Cutter):
	"""
	This class works in a similar way to L{Map <map.Map>}, except that it
	shows a scalar field cut using a plane. The plane can be translated 
	and rotated along the X, Y and Z axes.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):
		"""
		Initialise the MapOnPlanceCut.	

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
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme 
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		# Keeps track whether MapOnPlaneCut has been modified.
		self.__modified = True 
		Transform.__init__(self)	
		Plane.__init__(self)
		Cutter.__init__(self)
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

		# ----- Map on a plane -----
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
			c2p = CellDataToPointData(
					self.__data_collector._getDataCollectorOutput())
			self._setupCutter(c2p._getCellToPointOutput(), self._getPlane())
		elif(self.__cell_to_point == False): # No conversion happens.
			c2p = CellDataToPointData(
					self.__data_collector._getDataCollectorOutput())
			self._setupCutter(self.__data_collector._getDataCollectorOutput(),
				self._getPlane())

		self._setupDataSetMapper(self._getCutterOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the MapOnPlaneCut or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the surface map cut using a plane.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


###############################################################################


from clipper import Clipper

# NOTE: DataSetMapper, Actor3D, Transform, Plane and Clipper were inherited 
# to allow access to their public methods from the driver.
class MapOnPlaneClip(DataSetMapper, Actor3D, Transform, Plane, Clipper):
	"""
	This class works in a similar way to L{MapOnPlaneCut <map.MapOnPlaneCut>},
	except that it shows a scalar field clipped using a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):
		"""
		Initialise the MapOnPlaneClip.

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
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		# Keeps track whether MapOnPlaneClip has been modified.
		self.__modified = True 
		Transform.__init__(self)	
		Plane.__init__(self)
		Clipper.__init__(self)
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

		# ----- Map on a clipped plane -----

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
			c2p = CellDataToPointData(data_collector._getDataCollectorOutput())
			self._setupClipper(c2p._getCellToPointOutput(), self._getPlane())
		elif(self.__cell_to_point == False): # No conversion happens.
			self._setupClipper(data_collector._getDataCollectorOutput(),
					self._getPlane())

		self._setClipFunction()
		self._setupDataSetMapper(self._getClipperOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the MapOnPlaneClip or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the surface map clip using a plane.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


#############################################################################


# NOTE: DataSetMapper, Actor3D and Clipper were inherited 
# to allow access to their public methods from the driver.
class MapOnScalarClip(DataSetMapper, Actor3D, Clipper):
	"""
	This class works in a similar way to L{Map <map.Map>}, except that it
	shows a scalar field clipped using a scalar value.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):
		"""
		Initialise the MapOnScalarClip.

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
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		# Keeps track whether MapOnScalarClip has been modified.
		self.__modified = True 
		Clipper.__init__(self)
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

		# ----- Map clipped using a scalar value -----

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
			# None is used because a plane is not required when a scalar 
			# value is used to perform the clipping.
			self._setupClipper(c2p._getDataCollectorOutput(), None)
		elif(self.__cell_to_point == False): # No conversion happens.
			self._setupClipper(
					self.__data_collector._getDataCollectorOutput(),None)

		self._setupDataSetMapper(self._getClipperOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the MapOnScalarClip or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the surface map clip using scalar data.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False



