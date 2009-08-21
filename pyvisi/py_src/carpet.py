
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

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"

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

		:attention: The source can either be point or cell data. If the 		source is cell data, a conversion to point data may or may not be 		required, in order for the object to be rendered correctly. If a conversion is needed, the 'cell_to_point' flag must be set to 'True', otherwise 'False' (which is the default). On occasions, an inaccurate object may be rendered from cell data even after conversion. When 3D data is used, a cut will be performed on the scalar field using a plane before deformation occurs on the plane. However, if 2D data is used a cut will NOT be performed and deformation will instead occur immediately on the scalar field.  Pyvisi distinguishes 2D from 3D data by retrieving the length of the z-axis. A 2D data is assumed to have a z-axis length of zero while a 3D data is assumed to have a z-axis length of non-zero. There are exceptions to these rules where some 2D data may have a non-zero 		z-axis length. However, such exceptions are not taken into account 		at this stage.

		:type scene: `Scene` object
		:param scene: Scene in which objects are to be rendered on
		:type data_collector: `DataCollector` object
		:param data_collector: Deal with source of data for visualisation
		:type viewport: `Viewport` constant
		:param viewport: Viewport in which objects are to be rendered on
		:param warp_mode: `WarpMode` constant
		:type warp_mode: Mode in which to deform the scalar field 
		:type lut: `Lut` constant
		:param lut: Lookup table color scheme
		:type cell_to_point: Boolean
		:param cell_to_point: Converts cell data to point data (by averaging)
		:type outline: Boolean
		:param outline: Places an outline around the domain surface
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
			self.__lookup_table = LookupTable()
			self.__lookup_table._setTableValue()
		elif(self.__lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			self.__lookup_table = LookupTable()
			self.__lookup_table._setLookupTableToGreyScale()

		self._setupPlane(self._getTransform())

	def _isModified(self):	
		"""
		Return whether the Carpet or DataCollector has been modified.

		:rtype: Boolean
		:return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the carpet.

		:type scene: `Scene` object
		:param scene: Scene in which objects are to be rendered on
		"""

		# This entire 'if' section had to be moved from the __init__ method
		# due to the use of 'GetBounds()'. A source (i.e. xml file) must be
		# supplied before 'GetBounds()' is able to return the correct value
		# when executed. This is to accommodate for lazy evaluation.
		if(self.__modified):
			# Get the bounds of the object in the form of
			# (xmin, xmax, ymin, ymax, zmin, zmax).
			bounds = \
					self.__data_collector._getDataCollectorOutput().GetBounds()
			# Length of the z-axis (max - min). Assumption is made that if the
			# length of the z-axis is equal to zero, the the data set is 2D.
			# Otherwise, the data set is 3D. However, there are exceptions 
			# to this rule as some 2D data sets may have a z-axis length 
			# of non-zero, but such exceptions are not taken into account here.
			z_axis_length = bounds[5] - bounds[4]

			if(self.__cell_to_point == True): #Convert cell data to point data.
				c2p = CellDataToPointData(\
						self.__data_collector._getDataCollectorOutput())
				if(z_axis_length != 0): # A cutter is used for 3D data.
					self._setupCutter(c2p._getCellToPointOutput(), \
							self._getPlane()) 	
					self._setupWarp(self._getCutterOutput())
				elif(z_axis_length == 0): # A cutter is not used for 2D data.
					self._setupWarp(c2p._getCellToPointOutput())

			elif(self.__cell_to_point == False): # No conversion happens.	
				if(z_axis_length != 0): # A cutter is used for 3D data.
					self._setupCutter(\
							self.__data_collector._getDataCollectorOutput(), \
							self._getPlane()) 	
					self._setupWarp(self._getCutterOutput())
				elif(z_axis_length == 0): # A cutter is not used for 2D data.
					self._setupWarp(
							self.__data_collector._getDataCollectorOutput())

			self._setupDataSetMapper(self._getWarpOutput(),
					self.__lookup_table._getLookupTable())

			self._setupActor3D(self._getDataSetMapper())
			scene._addActor3D(self.__viewport, self._getActor3D())


		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()

			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(not(self._isScalarRangeSet())): 
				self._setScalarRange(self.__data_collector._getScalarRange())
			self.__modified = False



