#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

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
from actor import Actor2D, Actor3D
from lookuptable import LookupTable
from constant import Viewport, Color, Lut, LegendType
from scalarbar import ScalarBar

# NOTE: ScalarBarModule, DataSetMapper and Actor3D were inherited to allow 
# access to their public methods from the driver.
class Legend(ScalarBar, DataSetMapper, Actor3D):
	"""
	Class that shows a scalar field on a domain surface. The domain surface
	can either be colored or grey-scaled, depending on the lookup table used.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	def __init__(self, scene, data_collector, viewport = Viewport.SOUTH_WEST,\
			lut = Lut.COLOR, legend = LegendType.SCALAR):
		"""
		Initialise the Legend.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for vizualisation
		@type viewport: L{Viewport <constant.Viewport>} constant  
		@param viewport: Viewport in which objects are to be rendered on 
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type legend: L{Lut <constant.LegendType>} constant
		@param legend: Type of legend
		"""

		self.__data_collector = data_collector
		self.__viewport = viewport
		self.__lut = lut
		self.__legend = legend
		
		self.__modified = True # Keeps track whether Legend has been modified.
		ScalarBar.__init__(self)

		# NOTE: DataSetMapper and Actor3D were required in order for the 
		# scalar bar to be colored. If the mapper and actor were not used,
		# the scalar bar will be black in color.
		DataSetMapper.__init__(self)
		Actor3D.__init__(self)
		scene._addVisualizationModules(self)

		# ----- Scalar Bar -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(self.__lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable() 
			lookup_table._setTableValue()
		elif(self.__lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable() 
			lookup_table._setLookupTableToGreyScale()

		self._setupDataSetMapper(\
				self.__data_collector._getDataCollectorOutput(), \
				lookup_table._getLookupTable())	
		self._setScalarBarLookupTable(self._getDataSetMapperLookupTable())
		self.setOrientationToHorizontal()
		self._setupActor3D(self._getDataSetMapper())
		self.setOpacity(0)
		scene._addActor3D(self.__viewport, self._getActor3D())
		scene._addActor3D(self.__viewport, self._getScalarBar())
	
	def _isModified(self):	
		"""
		Return whether the Legend or DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__modified or self.__data_collector._isModified()

	def _render(self, scene):
		"""
		Render the legend.

		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		"""

		if (self._isModified() == True):
			if(self.__data_collector._isScalarSet() == True):
				self.__data_collector._setActiveScalar()

			# self._isScalarRangeSet checks whether the scalar range has been
			# specified by the user. If it has, then the scalar range
			# read from the source will be ignored.
			if(self.__legend == LegendType.SCALAR and \
					(not(self._isScalarRangeSet()))):
				self._setScalarRange(self.__data_collector._getScalarRange())
			elif(self.__legend == LegendType.VECTOR and \
					(not(self._isScalarRangeSet()))):
				self._setScalarRange(self.__data_collector._getVectorRange())
			self.__modified = False


