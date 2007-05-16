"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from actor import Actor3D
from lookuptable import LookupTable
from outline import Outline
from constant import Viewport, Color, Lut, VizType, ColorMode
from contourmodule import ContourModule
from average import CellDataToPointData

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

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the contour) can be executed. Actor3D 
		# methods associated with Outline cannot be invoked from the driver. 
		# They can only be called within here, which is why Outline must
		# be place before the contour as there is unlikely to be any changes
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

		# ----- Contour -----

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
			ContourModule.__init__(self, c2p._getOutput())	
		elif(cell_to_point == False): # No conversion happens.	
			ContourModule.__init__(self, data_collector._getOutput())	

		# By default 10 contours are generated and the scalar range is based
		# on the scalar data range.
		ContourModule.generateContours(self, 10, 
				data_collector._getScalarRange()[0],
				data_collector._getScalarRange()[1])

		DataSetMapper.__init__(self, ContourModule._getOutput(self), 
				lookup_table._getLookupTable())	
		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())

		data_collector._paramForUpdatingMultipleSources(VizType.CONTOUR,
				ColorMode.SCALAR, DataSetMapper._getDataSetMapper(self),
				ContourModule._getContour(self))

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))


###############################################################################


from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, ContourModule, Transform, Plane and Cutter were 
# inherited to allow access to their public methods from the driver.
class ContourOnPlaneCut(DataSetMapper, Actor3D, ContourModule, Transform, 
		Plane, Cutter):
	"""
	This class works in a similar way to L{MapOnPlaneCut <map.MapOnPlaneCut>},
	except that it shows a scalar field by contour surfaces on a plane.
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

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the contour) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must
		# be place before the contour as there is unlikely to be any changes
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

		# ----- Contour on a cut plane -----

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

		ContourModule.__init__(self, Cutter._getOutput(self))
		# By default 10 contours are generated and the scalar range is based
		# on the scalar data range.
		ContourModule.generateContours(self, 10, 
				data_collector._getScalarRange()[0],
				data_collector._getScalarRange()[1])

		DataSetMapper.__init__(self, ContourModule._getOutput(self), 
				lookup_table._getLookupTable())
		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())	

		data_collector._paramForUpdatingMultipleSources(VizType.CONTOUR,
				ColorMode.SCALAR, DataSetMapper._getDataSetMapper(self),
				ContourModule._getContour(self))

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))


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

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the contour) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must
		# be place before the contour as there is unlikely to be any changes
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

		# ----- Contour on a clipped plane -----

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
			Clipper.__init__(self, c2p._getOutput(), Plane._getPlane(self))
		elif(cell_to_point == False): # No conversion happens.	
			Clipper.__init__(self, data_collector._getOutput(), 
					Plane._getPlane(self))

		Clipper._setClipFunction(self)

		ContourModule.__init__(self, Clipper._getOutput(self))
		# By default 10 contours are generated and the scalar range is based
		# on the scalar data range.
		ContourModule.generateContours(self, 10, 
				data_collector._getScalarRange()[0],
				data_collector._getScalarRange()[1])

		DataSetMapper.__init__(self, ContourModule._getOutput(self), 
				lookup_table._getLookupTable())
		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())	

		data_collector._paramForUpdatingMultipleSources(VizType.CONTOUR,
				ColorMode.SCALAR, DataSetMapper._getDataSetMapper(self),
				ContourModule._getContour(self))

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))

