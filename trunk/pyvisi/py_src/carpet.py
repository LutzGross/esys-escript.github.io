"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut, WarpMode, VizType, ColorMode
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

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the carpet) can be executed. Actor3D 
		# methods associated with Outline cannot be invoked from the driver. 
		# They can only be called within here, which is why Outline must
		# be place before the carpet as there is unlikely to be any changes
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

		# ----- Carpet -----

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		Transform.__init__(self)
		Plane.__init__(self, Transform._getTransform(self))

		if(cell_to_point == True): # Converts cell data to point data.
			c2p = CellDataToPointData(data_collector._getOutput())
			Cutter.__init__(self, c2p._getOutput(), Plane._getPlane(self)) 	
		elif(cell_to_point == False): # No conversion happens.	
			Cutter.__init__(self, data_collector._getOutput(), 
					Plane._getPlane(self)) 	

		Warp.__init__(self, Cutter._getOutput(self), warp_mode)
			
		DataSetMapper.__init__(self, Warp._getOutput(self),
				lookup_table._getLookupTable())
		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())

		data_collector._paramForUpdatingMultipleSources(VizType.CARPET,
				ColorMode.SCALAR, DataSetMapper._getDataSetMapper(self))

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))

