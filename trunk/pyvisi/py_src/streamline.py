"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut, ColorMode, VizType
from streamlinemodule import  StreamLineModule
from tube import Tube
from point import PointSource
from outline import Outline
from average import CellDataToPointData

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

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Tube) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Tube as there is unlikely to be any changes
		# made to the Outline's Actor3D.

 		# ----- Outline -----

		if(outline == True):
			outline = Outline(data_collector._getOutput())
			DataSetMapper.__init__(self, outline._getOutput())

			Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
			# Default outline color is black.
			Actor3D.setColor(self, Color.BLACK)
			# Default line width is 1.
			Actor3D._setLineWidth(self, 1)
			scene._addActor3D(viewport, Actor3D._getActor3D(self))

		# ----- Streamline -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		if(cell_to_point == True): # Converts cell data to point data.
			c2p = CellDataToPointData(data_collector._getOutput())
			PointSource.__init__(self, c2p._getOutput())
			StreamLineModule.__init__(self, c2p._getOutput(), 
					PointSource._getOutput(self))
		elif(cell_to_point == False): # No conversion happens.
			PointSource.__init__(self, data_collector._getOutput())
			StreamLineModule.__init__(self, data_collector._getOutput(), 
					PointSource._getOutput(self))

		Tube.__init__(self, StreamLineModule._getOutput(self))
		DataSetMapper.__init__(self, Tube._getOutput(self), 
				lookup_table._getLookupTable())

		if(color_mode == ColorMode.VECTOR): # Color velocity by vector.
			DataSetMapper._setScalarVisibilityOn(self)
			StreamLineModule._setSpeedScalarsOn(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getVectorRange())

			data_collector._paramForUpdatingMultipleSources(VizType.STREAMLINE,
					ColorMode.VECTOR, DataSetMapper._getDataSetMapper(self))

		elif(color_mode == ColorMode.SCALAR): # Color velocity by scalar.
			DataSetMapper._setScalarVisibilityOn(self)
			StreamLineModule._setSpeedScalarsOff(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getScalarRange())
			
			data_collector._paramForUpdatingMultipleSources(VizType.STREAMLINE,
					ColorMode.SCALAR, DataSetMapper._getDataSetMapper(self))

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))

