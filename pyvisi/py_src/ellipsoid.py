"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut, ColorMode
from sphere import Sphere
from normals import Normals
from glyph import  TensorGlyph
from outline import Outline
from point import MaskPoints
from average import CellDataToPointData

# NOTE: DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph 
# and MaskPoints  were inherited to allow access to their 
# public methods from the driver.
class Ellipsoid(DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph, 
		MaskPoints):
	"""
	Class that shows a tensor field using ellipsoids. The ellipsoids can either
	be colored or grey-scaled, depending on the lookup table used.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True): 
		"""
		Initialise the Ellipsoid.

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

		# Keeps track whether Ellipsoid has been modified.
		self.__modified = True 		
		MaskPoints.__init__(self)
		Sphere.__init__(self)
		TensorGlyph.__init__(self) 
		Normals.__init__(self)
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

		# ----- Ellipsoid -----

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
			self._setupMaskPoints(
					self.__data_collector._getDataCollectorOutput())

		self._setupTensorGlyph(self._getMaskPointsOutput(), 
				self._getSphereOutput()) 
		self._setupNormals(self._getTensorGlyphOutput())

		self._setupDataSetMapper(self._getNormalsOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the Ellipsoid or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the ellipsoids.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			if(self.__data_collector._isTensorSet() == True):
				self.__data_collector._setActiveTensor()

			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


###############################################################################


from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph, Transform, Plane,
# Cutter and MaskPoints were inherited to allow access to 
# their public methods from the driver.
class EllipsoidOnPlaneCut(DataSetMapper, Actor3D, Sphere, Normals,  
		TensorGlyph, Transform, Plane, Cutter, MaskPoints):
	"""
	This class works in a similar way to L{MapOnPlaneCut <map.MapOnPlaneCut>},
	except that it shows a tensor field using ellipsoids cut using a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True): 
		"""
		Initialise the EllipsoidOnPlaneCut.

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

		# Keeps track whether EllipsoidOnPlaneCut has been modified.
		self.__modified = True 		
		Transform.__init__(self)
		Plane.__init__(self)
		Cutter.__init__(self)
		MaskPoints.__init__(self)
		Sphere.__init__(self)
		TensorGlyph.__init__(self) 
		Normals.__init__(self)
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

		# ----- Ellipsoid on a cut plane -----

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
			self._setupCutter(self.__data_collector._getDataCollectorOutput(), 
					self._getPlane()) 	

		self._setupMaskPoints(self._getCutterOutput())

		self._setupTensorGlyph(self._getMaskPointsOutput(),
				self._getSphereOutput())
		self._setupNormals(self._getTensorGlyphOutput()) 

		self._setupDataSetMapper(self._getNormalsOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the EllipsoidOnPlaneCut or DataCollector has been 
		modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the ellipsoids cut using a plane.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			if(self.__data_collector._isTensorSet() == True):
				self.__data_collector._setActiveTensor()
			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


###############################################################################


from clipper import Clipper

# NOTE: DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph, Transform, Plane,
# Clipper and MaskPoints were inherited to allow access to 
# their public methods from the driver.
class EllipsoidOnPlaneClip(DataSetMapper, Actor3D, Sphere, Normals,  
	TensorGlyph, Transform, Plane, Clipper, MaskPoints):
	"""
	This class works in a similar way to L{MapOnPlaneClip <map.MapOnPlaneClip>},
	except that it shows a tensor field using ellipsoids clipped using a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True): 
		"""
		Initialise the EllipsoidOnPlaneClip.

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
		@param viewport: Viewport in which object are to be rendered on
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

		# Keeps track whether EllipsoidOnPlaneClip has been modified.
		self.__modified = True 		
		Transform.__init__(self)
		Plane.__init__(self)
		Clipper.__init__(self)
		MaskPoints.__init__(self)
		Sphere.__init__(self)
		TensorGlyph.__init__(self) 
		Normals.__init__(self)
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

		# ----- Ellipsoid on a clipped plane -----

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

		# NOTE: TensorGlyph must come before Clipper. Otherwise clipping
		# may not work correctly.
		self._setupTensorGlyph(self._getMaskPointsOutput(),
				self._getSphereOutput())
		self._setupNormals(self._getTensorGlyphOutput())

		# NOTE: Clipper must come after TensorGlyph. Otherwise clipping
		# may not work correctly.
		self._setupClipper(self._getNormalsOutput(), 
				self._getPlane()) 	
		self._setClipFunction()

		self._setupDataSetMapper(self._getClipperOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())

	def _isModified(self):	
		"""
		Return whether the EllipsoidOnPlaneClip or DataCollector has been 
		modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the ellipsoids clip using a plane.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()
			if(self.__data_collector._isTensorSet() == True):
				self.__data_collector._setActiveTensor()
			self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


