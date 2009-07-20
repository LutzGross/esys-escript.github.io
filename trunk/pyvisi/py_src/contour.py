
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
from lookuptable import LookupTable
from outline import Outline
from constant import Viewport, Color, Lut, ColorMode
from contourmodule import ContourModule
from average import CellDataToPointData
from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

# NOTE: DataSetMapper, Actor3D and ContourModule were inherited to allow 
# access to their public methods from the driver.
class Contour(DataSetMapper, Actor3D, ContourModule):
	"""
	Class that shows a scalar field by contour surfaces. The contour surfaces
	can either be colored or grey-scaled, depending on the lookup table used.
	This class can also be used to generate iso surfaces.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):
		"""
		Initialise the Contour.

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
		
		self.__modified = True # Keeps track whether Contour has been modified.
		ContourModule.__init__(self)
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

		# ----- Contour -----

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
			self._setupContourModule(c2p._getCellToPointOutput())	
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupContourModule(
					self.__data_collector._getDataCollectorOutput())	

		self._setupDataSetMapper(self._getContourModuleOutput(), 
				lookup_table._getLookupTable())	

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())
	
	def _isModified(self):	
		"""
		Return whether the Contour or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the contour.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()

			# By default 10 contours are generated and the scalar range is 
			# based on the scalar data range.
			contours = 10
			lower_range = self.__data_collector._getScalarRange()[0]
			upper_range = self.__data_collector._getScalarRange()[1]

			if(self._isContoursSet() == True):
				contours = None	
			if(self._isLowerRangeSet() == True):
				lower_range = None
			if(self._isUpperRangeSet() == True):
				upper_range = None

			self.generateContours(contours, lower_range, upper_range)
			self._generateContours()
						
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

# NOTE: DataSetMapper, Actor3D, ContourModule, Transform, Plane and Cutter 
# were inherited to allow access to their public methods from the driver.
class ContourOnPlaneCut(DataSetMapper, Actor3D, ContourModule, Transform, 
		Plane, Cutter):
	"""
	This class works in a similar way to L{MapOnPlaneCut <map.MapOnPlaneCut>},
	except that it shows a scalar field by contour surfaces cut using a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):

		"""
		Initialise the ContourOnPlaneCut.

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
		
		# Keeps track whether ContourOnPlaneCut has been modified.
		self.__modified = True 		
		Transform.__init__(self)
		Plane.__init__(self)
		Cutter.__init__(self)
		ContourModule.__init__(self)
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

		# ----- Contour on a cut plane -----

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
			c2p = CellDataToPointData(self.__data_collector._getOutput())
			self._setupCutter(c2p._getCellToPointOutput(), self._getPlane())
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupCutter(self.__data_collector._getDataCollectorOutput(), 
					self._getPlane())

		self._setupContourModule(self._getCutterOutput())
		self._setupDataSetMapper(self._getContourModuleOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())
	
	def _isModified(self):	
		"""
		Return whether the ContourOnPlaneCut or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the contour cut using a plane.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()

			# By default 10 contours are generated and the scalar range is based
			# on the scalar data range.
			contours = 10
			lower_range = self.__data_collector._getScalarRange()[0]
			upper_range = self.__data_collector._getScalarRange()[1]

			if(self._isContoursSet() == True):
				contours = None	
			if(self._isLowerRangeSet() == True):
				lower_range = None
			if(self._isUpperRangeSet() == True):
				upper_range = None

			self.generateContours(contours, lower_range, upper_range)
			self._generateContours()
						
			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
				self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


###############################################################################


from clipper import Clipper

# NOTE: DataSetMapper, Actor3D, ContourModule, Transform, Plane and Clipper 
# were inherited to allow access to their public methods from the driver.
class ContourOnPlaneClip(DataSetMapper, Actor3D, ContourModule, Transform, 
		Plane, Clipper):
	"""
	This class works in a similar way to L{MapOnPlaneClip <map.MapOnPlaneClip>}
	, except that it shows a scalar field by contour surfaces clipped using
	a plane.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no lut is specified, the color scheme will be used. 
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST, 
			lut = Lut.COLOR, cell_to_point = False, outline = True):

		"""
		Initialise the ContourOnPlaneClip.

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
		
		# Keeps track whether ContourOnPlaneClip has been modified.
		self.__modified = True 		
		Transform.__init__(self)
		Plane.__init__(self)
		Clipper.__init__(self)
		ContourModule.__init__(self)
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

		# ----- Contour on a clipped plane -----

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
			self._setupClipper(c2p._getCellToPointOutput(), self._getPlane())
		elif(self.__cell_to_point == False): # No conversion happens.	
			self._setupClipper(self.__data_collector._getDataCollectorOutput(), 
					self._getPlane())

		self._setClipFunction()
		self._setupContourModule(self._getClipperOutput())

		self._setupDataSetMapper(self._getContourModuleOutput(), 
				lookup_table._getLookupTable())

		self._setupActor3D(self._getDataSetMapper())
		scene._addActor3D(self.__viewport, self._getActor3D())
	
	def _isModified(self):	
		"""
		Return whether the ContourOnPlaneClip or DataCollector has been 
		modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the contour clip using a plane.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()

			# By default 10 contours are generated and the scalar range is based
			# on the scalar data range.
			contours = 10
			lower_range = self.__data_collector._getScalarRange()[0]
			upper_range = self.__data_collector._getScalarRange()[1]

			if(self._isContoursSet() == True):
				contours = None	
			if(self._isLowerRangeSet() == True):
				lower_range = None
			if(self._isUpperRangeSet() == True):
				upper_range = None

			self.generateContours(contours, lower_range, upper_range)
			self._generateContours()
						
			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
				self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False


