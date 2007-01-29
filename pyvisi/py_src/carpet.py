"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut
from warp import Warp
from outline import Outline
from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, Warp, Transform, Plane and Cutter  were 
# inherited to allow access to their public methods from the driver.
class Carpet(DataSetMapper, Actor3D, Warp, Transform, 
		Plane, Cutter):
	"""
	Class that shows a scalar field on a plance deformated along the place 
	normal.
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no scalar field is specified, the first encountered in the file will
	# be loaded automatically. If no lut is specified, the color scheme will
	# be used. 
	def __init__(self, scene, data_collector, scalar = None, warp_mode = None, 
			viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type scalar: String
		@param scalar: Scalar field to load from the source file
		@param warp_mode: L{WarpMode {<constant.WarpMode>} constant
		@type warp_mode: Mode in which to deform the scalar data
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

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

		if(scalar != None):
			data_collector._setActiveScalar(scalar)
				
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

		Cutter.__init__(self, data_collector._getOutput(),
				Plane._getPlane(self))

		#warp = Warp(Cutter._getOutput(self), warp_mode)
		Warp.__init__(self, Cutter._getOutput(self), warp_mode)
			
		#DataSetMapper.__init__(self, warp._getOutput(),
		DataSetMapper.__init__(self, Warp._getOutput(self),
				lookup_table._getLookupTable())

		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))

