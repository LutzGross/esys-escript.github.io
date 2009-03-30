
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
from actor import Actor3D
from lookuptable import LookupTable
from outline import Outline
from constant import Viewport, Color, Lut, ColorMode
from average import CellDataToPointData
from esys.escript import getMPISizeWorld

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
		'True', otherwise 'False' (which is the default). On occasions, an
		inaccurate object may be rendered from cell data even after conversion.

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

			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
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
		'True', otherwise 'False' (which is the default). On occasions, an
		inaccurate object may be rendered from cell data even after conversion.

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
			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
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
		'True', otherwise 'False' (which is the default). On occasions, an
		inaccurate object may be rendered from cell data even after conversion.

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
			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
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
		'True', otherwise 'False' (which is the default). On occasions, an
		inaccurate object may be rendered from cell data even after conversion.

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
			self._setupClipper(c2p._getCellToPointOutput(), None)
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
			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
				self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


##############################################################################


from rotation import Rotation
from geometry import Geometry

# NOTE: DataSetMapper, Actor3D, Clipper and Rotation were inherited 
# to allow access to their public methods from the driver.
class MapOnScalarClipWithRotation(DataSetMapper, Actor3D, Clipper, Rotation):
	"""
	This class works in a similar way to L{Map <map.Map>}, except that it
	shows a 2D scalar field clipped using a scalar value and subsequently
	rotated around the z-axis to create a  3D looking effect. This class should 
	only be used with 2D data sets and NOT 3D.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False):
		"""
		Initialise the MapOnScalarClipWithRotation.

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
		@param viewport: Viewport in which objects are to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		
		# Keeps track whether MapOnScalarClipWithRotation has been modified.
		self.__modified = True 
		Clipper.__init__(self)
		Rotation.__init__(self)
		DataSetMapper.__init__(self)
		Actor3D.__init__(self)
		scene._addVisualizationModules(self)

	def _isModified(self):	
		"""
		Return whether the MapOnScalarClipWithRotation or DataCollector has 
		been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the surface map clip using scalar data and subsequently rotated.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		# This entire 'if' section had to be moved from the __init__ method
		# due to the use of 'GetPoints().GetNumberOfPoints()'. A source 
		# (i.e. xml file) must be 
		# supplied before 'GetPoints().GetNumberOfPoints()' is able to return 
		# the correct value
		# when executed. This is to accommodate for lazy evaluation.
		if (self._isModified() == True):

			# Swaps the points from the y-axis to the z-axis and vice-versa.
			# This is needed as VTK is only able to perform rotation along
			# the z-axis.
			output = self.__data_collector._getDataCollectorOutput()
			num_points = output.GetPoints().GetNumberOfPoints()
			for i in range (num_points):
				point = output.GetPoints().GetPoint(i)
				output.GetPoints().SetPoint(i, point[0], point[2], point[1])

			# --- Map clipped using a scalar value and subsequently rotated ---

			# NOTE: Lookup table color mapping (color or grey scale) MUST be set
			# before DataSetMapper. If it is done after DataSetMapper, no effect
			# will take place.
			if(self.__lut == Lut.COLOR): # Colored lookup table.
				lookup_table = LookupTable() 
				lookup_table._setTableValue()
			elif(self.__lut == Lut.GREY_SCALE): # Grey scaled lookup table.
				lookup_table = LookupTable() 
				lookup_table._setLookupTableToGreyScale()

			if(self.__cell_to_point == True): #Converts cell data to point data.
				c2p = CellDataToPointData(output)
				# None is used because a plane is not required when a scalar 
				# value is used to perform the clipping.
				self._setupClipper(c2p._getCellToPointOutput(), None)
			elif(self.__cell_to_point == False): # No conversion happens.
				self._setupClipper(output, None)

			geometry = Geometry(self._getClipperOutput())	
			
			self._setupRotationExtrusionFilter(geometry._getGeometryOutput())
			self._setupDataSetMapper(self._getRotationOutput(),
					lookup_table._getLookupTable())

			self._setupActor3D(self._getDataSetMapper())
			scene._addActor3D(self.__viewport, self._getActor3D())

			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
				self._setScalarRange(self.__data_collector._getScalarRange())

			self.__modified = False



