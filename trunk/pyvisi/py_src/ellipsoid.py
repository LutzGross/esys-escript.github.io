"""
@author: John NGUI
"""

import vtk
from mapper import DataSetMapper
from lookuptable import LookupTable
from actor import Actor3D
from constant import Viewport, Color, Lut
from sphere import Sphere
from normals import Normals
from glyph import  TensorGlyph
from outline import Outline

# NOTE: DataSetMapper, Actor3D, Sphere, Normals and TensorGlyph were inherited 
# to allow access to their public methods from the driver.
class Ellipsoid(DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph):
	"""
	Class that show a tensor field using ellipsoid.	
	"""
	
	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no vector field is specified, the first encountered in the file will
	# be loaded automatically.
	def __init__(self, scene, data_collector, tensor = None, 
			viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True): 

		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type tensor: String
		@param tensor: Tensor field to load from the source file
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Ellipsoid) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Ellipsoid as there is unlikely to be any changes
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


		# ----- Ellipsoid -----

		if(tensor != None):
			data_collector._setActiveTensor(tensor)

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		Sphere.__init__(self)
		TensorGlyph.__init__(self, data_collector._getOutput(), 
				Sphere._getOutput(self)) 
		Normals.__init__(self, TensorGlyph._getOutput(self))

		DataSetMapper.__init__(self, Normals._getOutput(self), 
				lookup_table._getLookupTable())

		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))


from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph, Transform, Plane,
# and Cutter were inherited to allow access to their public methods from 
# the driver.
class EllipsoidOnPlaneCut(DataSetMapper, Actor3D, Sphere, Normals,  
	TensorGlyph, Transform, Plane, Cutter):
	"""
	Class that show a tensor field using ellipsoids on a plane.	
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no vector field is specified, the first encountered in the file will
	# be loaded automatically.
	def __init__(self, scene, data_collector, tensor = None, 
			viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True): 

		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type tensor: String
		@param tensor: Tensor field to load from the source file
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Ellipsoid) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Ellipsoid as there is unlikely to be any changes
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


		# ----- Ellipsoid on a plane -----

		if(tensor != None):
			data_collector._setActiveTensor(tensor)

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

		Sphere.__init__(self)
		TensorGlyph.__init__(self, Cutter._getOutput(self),
				Sphere._getOutput(self))
		Normals.__init__(self, TensorGlyph._getOutput(self))

		DataSetMapper.__init__(self, Normals._getOutput(self), 
			lookup_table._getLookupTable())

		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())
		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))


from clipper import Clipper

# NOTE: DataSetMapper, Actor3D, Sphere, Normals, TensorGlyph, Transform, Plane,
# and Clipper were inherited to allow access to their public methods from 
# the driver.
class EllipsoidOnPlaneClip(DataSetMapper, Actor3D, Sphere, Normals,  
	TensorGlyph, Transform, Plane, Clipper):
	"""
	Class that show a tensor field using ellipsoids on a plane.	
	"""

	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no vector field is specified, the first encountered in the file will
	# be loaded automatically.
	def __init__(self, scene, data_collector, tensor = None, 
			viewport = Viewport.SOUTH_WEST, lut = Lut.COLOR, outline = True): 

		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type tensor: String
		@param tensor: Tensor field to load from the source file
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		"""

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Ellipsoid) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Ellipsoid as there is unlikely to be any changes
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


		# ----- Ellipsoid on a clipped plane -----

		if(tensor != None):
			data_collector._setActiveTensor(tensor)

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

		# NOTE: TensorGlyph must come before Clipper. Otherwise the output
		# will be incorrect.
		Sphere.__init__(self)
		TensorGlyph.__init__(self, data_collector._getOutput(),
				Sphere._getOutput(self))
		Normals.__init__(self, TensorGlyph._getOutput(self))


		# NOTE: Clipper must come after TensorGlyph. Otherwise the output 
		# will be incorrect.
		Clipper.__init__(self, Normals._getOutput(self), 
				Plane._getPlane(self)) 	
		Clipper._setClipFunction(self)


		DataSetMapper.__init__(self, Clipper._getOutput(self), 
			lookup_table._getLookupTable())

		DataSetMapper._setScalarRange(self, data_collector._getScalarRange())
		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))

