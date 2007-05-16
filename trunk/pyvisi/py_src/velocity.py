"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Arrow, ColorMode, Lut 
from arrow import Arrow2D, Arrow3D
from glyph import  Glyph3D
from outline import Outline
from point import MaskPoints
from average import CellDataToPointData

# NOTE: DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D and
# MaskPoints were inherited to allow access to their public 
# methods from the driver.
class Velocity(DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D, MaskPoints):
	"""
	Class that shows a vector field using arrows. The arrows can either be
	colored or grey-scaled, depending on the lookup table used. If the arrows
	are colored, there are two possible coloring modes, either using vector data
	or scalar data. Similarly, there are two possible types of
	arrows, either using two-dimensional or three-dimensional.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, arrow = Arrow.TWO_D,
			color_mode = ColorMode.VECTOR, viewport = Viewport.SOUTH_WEST,  
			lut = Lut.COLOR, cell_to_point = False, outline = True): 
		"""
		Initialise the Velocity.

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
		@type arrow: L{Arrow <constant.Arrow>} constant 
		@param arrow: Type of arrow (two dimensional or three dimensional)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__scene = scene
		self.__data_collector = data_collector
		self.__arrow = arrow
		self.__color_mode = color_mode
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		self.__modified = True # Keeps track whether Velocity has been modified.
		MaskPoints.__init__(self)
		Glyph3D.__init__(self) 
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
			self.__scene._addActor3D(self.__viewport, actor3D._getActor3D())

		# ----- Velocity -----

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
			self._setupMaskPoints(c2p._getCellToPointOutput())
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupMaskPoints(data_collector._getDataCollectorOutput())
	
		if(self.__arrow == Arrow.TWO_D): # Use 2D arrows.
			Arrow2D.__init__(self)
			self._setupGlyph3D(self._getMaskPointsOutput(), 
					self._getArrow2DOutput()) 
		elif(arrow == Arrow.THREE_D): # Use 3D arrows.
			Arrow3D.__init__(self)
			self._setupGlyph3D(self._getMaskPointsOutput(), 
					self._getArrow3DOutput()) 

		self._setupDataSetMapper(self._getGlyph3DOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		self.__scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the Velocity or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self):
		"""
		Render the velocity.
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			if(self.__data_collector._isVectorSet() == True):
				self.__data_collector._setActiveVector()

			# Color velocity by vector.
			if(self.__color_mode == ColorMode.VECTOR): 				
				self._setColorModeByVector()
				self._setRange(self.__data_collector._getVectorRange())
				self._setScalarRange(self.__data_collector._getVectorRange())
			# Color velocity by scalar.
			elif(self.__color_mode == ColorMode.SCALAR): 				
				self._setColorModeByScalar()
				self._setRange(self.__data_collector._getScalarRange())
				self._setScalarRange(self.__data_collector._getScalarRange())

			self.__modified = False


###############################################################################


from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D, Transform, Plane,
# Cutter and MaskPoints were inherited to allow access to 
# their public methods from the driver.
class VelocityOnPlaneCut(DataSetMapper, Actor3D, Arrow2D, Arrow3D,  
		Glyph3D, Transform, Plane, Cutter, MaskPoints):
	"""
	This class works in a similar way to L{MapOnPlaneCut <map.MapOnPlaneCut>},
	except that it shows a vector field using arrows cut using a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme willbe used.
	def __init__(self, scene, data_collector, arrow = Arrow.TWO_D, 
			color_mode = ColorMode.VECTOR, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True): 
		"""
		Initialise the VelocityOnPlaneCut.

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
		@type arrow: L{Arrow <constant.Arrow>} constant 
		@param arrow: Type of arrow (two dimensional or three dimensional)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__scene = scene
		self.__data_collector = data_collector
		self.__arrow = arrow
		self.__color_mode = color_mode
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		# Keeps track whether VelocityOnPlaneCut has been modified.
		self.__modified = True 
		Transform.__init__(self)
		Plane.__init__(self)
		Cutter.__init__(self)
		MaskPoints.__init__(self)
		Glyph3D.__init__(self) 
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
			self.__scene._addActor3D(self.__viewport, actor3D._getActor3D())

		# ----- Velocity on a cut plane -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		self._setupPlane(self._getTransform())

		if(self.__cell_to_point == True): # Converts cell data to point data.
			c2p = CellDataToPointData(
					self.__data_collector._getDataCollectorOutput())
			self._setupCutter(c2p._getCellToPointOutput(), self._getPlane()) 	
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupCutter(self.__data_collector._getDataCollectorOutput(), 
					self._getPlane()) 	

		self._setupMaskPoints(self._getCutterOutput())

		if(self.__arrow == Arrow.TWO_D): # Use 2D arrows.
			Arrow2D.__init__(self)
			self._setupGlyph3D(self._getMaskPointsOutput(), 
					self._getActor2DOutput()) 
		elif(self.__arrow == Arrow.THREE_D): # Use 3D arrows.
			Arrow3D.__init__(self)
			self._setupGlyph3D(self._getMaskPointsOutput(), 
					self._getArrow3DOutput()) 

		self._setupDataSetMapper(self._getGlyph3DOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		self.__scene._addActor3D(self.__viewport, self._getActor3D())
	
	def _isModified(self):	
		"""
		Return whether the VelocityOnPlaneCut or DataCollector has been 
		modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self):
		"""
		Render the velocity cut using a plane..
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isVectorSet() == True):
				self.__data_collector._setActiveVector()
			# Color velocity by vector.
			if(self.__color_mode == ColorMode.VECTOR): 				
				self._setColorModeByVector()
				self._setRange(self.__data_collector._getVectorRange())
				self._setScalarRange(self.__data_collector._getVectorRange())
			# Color velocity by scalar.
			elif(self.__color_mode == ColorMode.SCALAR): 				
				self._setColorModeByScalar()
				self._setRange(self.__data_collector._getScalarRange())
				self._setScalarRange(self.__data_collector._getScalarRange())

			self.__modified = False


###############################################################################


from clipper import Clipper

# NOTE: DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D, Transform, Plane,
# Clipper and MaskPoints  were inherited to allow access to 
# their public methods from the driver.
class VelocityOnPlaneClip(DataSetMapper, Actor3D, Arrow2D, Arrow3D,  
		Glyph3D, Transform, Plane, Clipper, MaskPoints):
	"""
	This class works in a similar way to L{MapOnPlaneClip <map.MapOnPlaneClip>}
	, except that it shows a vector field using arrows clipped using a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used.
	def __init__(self, scene, data_collector, arrow = Arrow.TWO_D, 
			color_mode = ColorMode.VECTOR, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True): 
		"""
		Initialise the VelocityOnPlaneClip.

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
		@type arrow: L{Arrow <constant.Arrow>} constant 
		@param arrow: Type of arrow (two dimensional or three dimensional)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which objects are to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type cell_to_point: Boolean
		@param cell_to_point: Converts cell data to point data (by averaging)
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		self.__scene = scene
		self.__data_collector = data_collector
		self.__arrow = arrow
		self.__color_mode = color_mode
		self.__viewport = viewport
		self.__lut = lut
		self.__cell_to_point = cell_to_point
		self.__outline = outline
		
		# Keeps track whether VelocityOnPlaneCut has been modified.
		self.__modified = True 
		Transform.__init__(self)
		Plane.__init__(self)
		Clipper.__init__(self)
		MaskPoints.__init__(self)
		Glyph3D.__init__(self) 
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
			self.__scene._addActor3D(self.__viewport, actor3D._getActor3D())

		# ----- Velocity on a clipped plane -----

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
			self._setupMaskPoints(c2p._getCellToPointOutput())
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupMaskPoints(
					self.__data_collector._getDataCollectorOutput())

		# NOTE: Glyph3D must come before Clipper. Otherwise clipping may not
		# work correctly.
		if(self.__arrow == Arrow.TWO_D): # Use 2D arrows.
			Arrow2D.__init__(self)
			self._setupGlyph3D(self._getMaskPointsOutput(), 
					self._getArrow2DOutput()) 
		elif(self.__arrow == Arrow.THREE_D): # Use 3D arrows.
			Arrow3D.__init__(self)
			self._setupGlyph3D(self._getMaskPointsOutput(), 
					self._getArrow3DOutput()) 
		
		self._setupClipper(self._getGlyph3DOutput(), self._getPlane()) 	
		self._setClipFunction()

		# NOTE: Clipper must come after Glyph. Otherwise clipping
        # may not work correctly.
		self._setupDataSetMapper(self._getClipperOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		self.__scene._addActor3D(self.__viewport, self._getActor3D())
	
	def _isModified(self):	
		"""
		Return whether the VelocityOnPlaneClip or DataCollector has been 
		modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self):
		"""
		Render the velocity clip using a plane..
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isVectorSet() == True):
				self.__data_collector._setActiveVector()
			# Color velocity by vector.
			if(self.__color_mode == ColorMode.VECTOR): 				
				self._setColorModeByVector()
				self._setRange(self.__data_collector._getVectorRange())
				self._setScalarRange(self.__data_collector._getVectorRange())
			# Color velocity by scalar.
			elif(self.__color_mode == ColorMode.SCALAR): 				
				self._setColorModeByScalar()
				self._setRange(self.__data_collector._getScalarRange())
				self._setScalarRange(self.__data_collector._getScalarRange())

			self.__modified = False


