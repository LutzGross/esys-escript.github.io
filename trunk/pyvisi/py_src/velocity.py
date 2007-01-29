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
from probe import Probe
from point import StructuredPoints 

# NOTE: DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D, StructuredPoints and
# Probe were inherited to allow access to their public methods from the driver.
class Velocity(DataSetMapper, Actor3D, Arrow2D, Arrow3D,  Glyph3D, StructuredPoints, Probe):
	"""
	Class that show a vector field using arrows.	
	"""
	
	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no vector field is specified, the first encountered in the file will
	# be loaded automatically.
	def __init__(self, scene, data_collector, vector = None, 
			viewport = Viewport.SOUTH_WEST, outline = True, 
			arrow = Arrow.TWO_D, color_mode = ColorMode.VECTOR, 
			lut = Lut.COLOR): 

		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type vector: String
		@param vector: Vector field to load from the source file
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		@type arrow: L{Arrow <constant.Arrow>} constant 
		@param arrow: Type of arrow (two dimensional or three dimensional)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		"""

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Velocity) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Velocity as there is unlikely to be any changes
		# made to the Outline's Actor3D.

 		# ----- Outline -----

		if(outline == True):
			#outline = Outline(Glyph3D._getOutput(self))
			outline = Outline(data_collector._getOutput())
			DataSetMapper.__init__(self, outline._getOutput())

			Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
			# Default outline color is black.
			Actor3D.setColor(self, Color.BLACK)
			# Default line width is 1.
			Actor3D._setLineWidth(self, 1)
			scene._addActor3D(viewport, Actor3D._getActor3D(self))


		# ----- Velocity -----

		if(vector != None):
			data_collector._setActiveVector(vector)

		# NOTE: Lookup table color mapping (color or grey scale) MUST be set
		# before DataSetMapper. If it is done after DataSetMapper, no effect
		# will take place.
		if(lut == Lut.COLOR): # Colored lookup table.
			lookup_table = LookupTable()
			lookup_table._setTableValue()
		elif(lut == Lut.GREY_SCALE): # Grey scaled lookup table.
			lookup_table = LookupTable()
			lookup_table._setLookupTableToGreyScale()

		StructuredPoints.__init__(self, data_collector._getOutput())
		Probe.__init__(self, data_collector._getOutput(), 
				StructuredPoints._getStructuredPoints(self))


		if(arrow == Arrow.TWO_D): # Use 2D arrows.
			Arrow2D.__init__(self)
			#Glyph3D.__init__(self, data_collector._getOutput(), 
			Glyph3D.__init__(self, Probe._getOutput(self), 
			#Glyph3D.__init__(self, MaskPoints._getOutput(self), 
					Arrow2D._getOutput(self), data_collector._getScalarRange()) 
		elif(arrow == Arrow.THREE_D): # Use 3D arrows.
			Arrow3D.__init__(self)
			#Glyph3D.__init__(self, data_collector._getOutput(), 
			Glyph3D.__init__(self, Probe._getOutput(self), 
			#Glyph3D.__init__(self, MaskPoints._getOutput(self), 
					Arrow3D._getOutput(self), data_collector._getScalarRange()) 

		DataSetMapper.__init__(self, Glyph3D._getOutput(self), 
				lookup_table._getLookupTable())

		if(color_mode == ColorMode.VECTOR): # Color velocity by vector.
			Glyph3D._setColorModeByVector(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getVectorRange())
		elif(color_mode == ColorMode.SCALAR): # Color velocity by scalar.
			Glyph3D._setColorModeByScalar(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getScalarRange())

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))


from transform import Transform
from plane import Plane
from cutter import Cutter

# NOTE: DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D, Transform, Plane,
# Cutter, StructuredPoints and  Probe were inherited to allow access to 
# their public methods from the driver.
class VelocityOnPlaneCut(DataSetMapper, Actor3D, Arrow2D, Arrow3D,  
		Glyph3D, Transform, Plane, Cutter, StructuredPoints, Probe):
	"""
	Class that show a vector field using arrows on a plane.	
	"""
	
	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no vector field is specified, the first encountered in the file will
	# be loaded automatically.
	def __init__(self, scene, data_collector, vector = None, 
			viewport = Viewport.SOUTH_WEST, outline = True, 
			arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR,
			lut = Lut.COLOR): 

		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type vector: String
		@param vector: Vector field to load from the source file
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		@type arrow: L{Arrow <constant.Arrow>} constant 
		@param arrow: Type of arrow (two dimensional or three dimensional)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		"""

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Velocity) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Velocity as there is unlikely to be any changes
		# made to the Outline's Actor3D.

 		# ----- Outline -----

		if(outline == True):
			#outline = Outline(Glyph3D._getOutput(self))
			outline = Outline(data_collector._getOutput())
			DataSetMapper.__init__(self, outline._getOutput())

			Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
			# Default outline color is black.
			Actor3D.setColor(self, Color.BLACK)
			# Default line width is 1.
			Actor3D._setLineWidth(self, 1)
			scene._addActor3D(viewport, Actor3D._getActor3D(self))


		# ----- Velocity on a cut plane -----

		if(vector != None):
			data_collector._setActiveVector(vector)

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

		StructuredPoints.__init__(self, data_collector._getOutput())
		Probe.__init__(self, data_collector._getOutput(), 
				StructuredPoints._getStructuredPoints(self))

		Cutter.__init__(self, Probe._getOutput(self), 
		#Cutter.__init__(self, data_collector._getOutput(), 
				Plane._getPlane(self)) 	

		if(arrow == Arrow.TWO_D): # Use 2D arrows.
			Arrow2D.__init__(self)
			Glyph3D.__init__(self, Cutter._getOutput(self), 
					#Arrow2D._getOutput(self), data_collector._getVectorRange()) 
					Arrow2D._getOutput(self), data_collector._getScalarRange()) 
		elif(arrow == Arrow.THREE_D): # Use 3D arrows.
			Arrow3D.__init__(self)
			Glyph3D.__init__(self, Cutter._getOutput(self), 
					#Arrow3D._getOutput(self), data_collector._getVectorRange()) 
					Arrow3D._getOutput(self), data_collector._getScalarRange()) 

		DataSetMapper.__init__(self, Glyph3D._getOutput(self), 
				lookup_table._getLookupTable())

		if(color_mode == ColorMode.VECTOR): # Color velocity by vector.
			Glyph3D._setColorModeByVector(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getVectorRange())
		elif(color_mode == ColorMode.SCALAR): # Color velocity by scalar.
			Glyph3D._setColorModeByScalar(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getScalarRange())

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))


from clipper import Clipper

# NOTE: DataSetMapper, Actor3D, Arrow2D, Arrow3D, Glyph3D, Transform, Plane,
# Clipper, StructuredPoints and Probe  were inherited to allow access to 
# their public methods from the driver.
class VelocityOnPlaneClip(DataSetMapper, Actor3D, Arrow2D, Arrow3D,  
		Glyph3D, Transform, Plane, Clipper, StructuredPoints, Probe):
	"""
	Class that show a vector field using arrows on a clipped plane.	
	"""
	
	# The SOUTH_WEST default viewport is used when there is only one viewport.
	# This saves the user from specifying the viewport when there is only one.
	# If no vector field is specified, the first encountered in the file will
	# be loaded automatically.
	def __init__(self, scene, data_collector, vector = None, 
			viewport = Viewport.SOUTH_WEST, outline = True, 
			arrow = Arrow.THREE_D, color_mode = ColorMode.VECTOR, 
			lut = Lut.COLOR): 

		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which objects are to be rendered on
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
				object
		@param data_collector: Deal with source of data for visualisation
		@type vector: String
		@param vector: Vector field to load from the source file
		@type viewport: L{Viewport <constant.Viewport>} constant
		@param viewport: Viewport in which the object is to be rendered on
		@type outline: Boolean
		@param outline: Places an outline around the domain surface
		@type arrow: L{Arrow <constant.Arrow>} constant 
		@param arrow: Type of arrow (two dimensional or three dimensional)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type lut : L{Lut <constant.Lut>} constant
		@param lut: Lookup table color scheme
		"""

		# NOTE: Actor3D is inherited and there are two instances declared here.
		# As a result, when methods from Actor3D is invoked from the driver,
		# only the methods associated with the latest instance (which in this
		# case is the Actor3D for the Velocity) can be executed. Actor3D
		# methods associated with Outline cannot be invoked from the driver.
		# They can only be called within here, which is why Outline must be
		# place before Velocity as there is unlikely to be any changes
		# made to the Outline's Actor3D.

 		# ----- Outline -----

		if(outline == True):
			#outline = Outline(Glyph3D._getOutput(self))
			outline = Outline(data_collector._getOutput())
			DataSetMapper.__init__(self, outline._getOutput())

			Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
			# Default outline color is black.
			Actor3D.setColor(self, Color.BLACK)
			# Default line width is 1.
			Actor3D._setLineWidth(self, 1)
			scene._addActor3D(viewport, Actor3D._getActor3D(self))


		# ----- Velocity on a clipped plane -----

		if(vector != None):
			data_collector._setActiveVector(vector)

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

		StructuredPoints.__init__(self, data_collector._getOutput())
		Probe.__init__(self, data_collector._getOutput(), 
				StructuredPoints._getStructuredPoints(self))

		#MaskPoints.__init__(self, data_collector._getOutput())

		# NOTE: Glyph3D must come before Clipper. Otherwise, the output will
		# be incorrect.
		if(arrow == Arrow.TWO_D): # Use 2D arrows.
			Arrow2D.__init__(self)
			#Glyph3D.__init__(self, data_collector._getOutput(), 
			Glyph3D.__init__(self, Probe._getOutput(self), 
			#Glyph3D.__init__(self, MaskPoints._getOutput(self), 
					Arrow2D._getOutput(self), data_collector._getScalarRange()) 
		elif(arrow == Arrow.THREE_D): # Use 3D arrows.
			Arrow3D.__init__(self)
			#Glyph3D.__init__(self, data_collector._getOutput(), 
			Glyph3D.__init__(self, Probe._getOutput(self), 
			#Glyph3D.__init__(self, MaskPoints._getOutput(self), 
					Arrow3D._getOutput(self), data_collector._getScalarRange()) 

		
		# NOTE: Clipper must come after Glyph3D. Otherwise, the output will
		# be incorrect.
		Clipper.__init__(self, Glyph3D._getOutput(self), 
				Plane._getPlane(self)) 	
		Clipper._setClipFunction(self)

		DataSetMapper.__init__(self, Clipper._getOutput(self), 
				lookup_table._getLookupTable())

		if(color_mode == ColorMode.VECTOR): # Color velocity by vector.
			Glyph3D._setColorModeByVector(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getVectorRange())
		elif(color_mode == ColorMode.SCALAR): # Color velocity by scalar.
			Glyph3D._setColorModeByScalar(self)
			DataSetMapper._setScalarRange(self, 
					data_collector._getScalarRange())

		Actor3D.__init__(self, DataSetMapper._getDataSetMapper(self))
		scene._addActor3D(viewport, Actor3D._getActor3D(self))
