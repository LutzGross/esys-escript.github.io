"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common 
from geo import Position

class Plane(Common):
	"""
	Class that performs cutting and clipping using a plane 
	as its implicit function.
	"""

	def __init__(self, scene, data_collector, component, transform, lut, mode): 
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@type component: String
		@param component: Component to be cut using the plane
		@type transform: L{Transform <geo.Transform>} object
		@param transform: Orientation of the plane
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		Common.__init__(self, scene, data_collector)
		if(mode == "Cut" or mode == "Clip"):
			self.vtk_plane = vtk.vtkPlane() 
			self.transform = transform.getTransform() 
			self.setPlane()

		# Convert unstructured grid data to polygonal data.
		#vtk_geometry = vtk.vtkGeometryFilter()
		#vtk_geometry.SetInput(component)
	
		if(mode == "Cut"): # Executed only when cutting is performed.
			self.vtk_cutter = vtk.vtkCutter()
			self.setCutter(component)
			Common.setMapperInput(self, self.vtk_cutter.GetOutput(), lut)
		elif(mode == "Clip"): # Executed only when clipping is performed.
			#self.vtk_clipper = vtk.vtkClipPolyData()
			self.vtk_clipper = vtk.vtkClipDataSet()
			#self.setClipper(vtk_geometry.GetOutput())
			self.setClipper(component)
			Common.setMapperInput(self, self.vtk_clipper.GetOutput(), lut)
		elif(mode == "ScalarClip"): # Executed only for scalar clipping 
			self.vtk_clipper = vtk.vtkClipDataSet()
			self.vtk_clipper.SetInput(component)
			self.setInsideOutOn()
			Common.setMapperInput(self, self.vtk_clipper.GetOutput(), lut)

		Common.setActorInput(self)
		Common.addActor(self)

	def setPlane(self):
		"""
		Set up the plane.
		"""

		#print "--- in set plane..."	
		#bound = self.data_collector.getReader().GetOutput().GetBounds()
		#print bound
		#print bound[0]
		#print bound[1]
		#print bound[2]
		#print bound[3]
		#print bound[4]
		#print bound[5]
	
		
		# Default plane is an XY plane
		# Default plane origin is the center. This causes the plane to cut 
		# and clip at the center to the rendered object.
		#self.vtk_plane.SetOrigin(
		#	self.data_collector.getReader().GetOutput().GetCenter())
		# Default plane origin. This causes the plane to cut and clip at 
		# the base/origin of the rendered object.
		self.setPlaneOrigin(Position(0,0,0))
		# Default plane normal is orthorgonal to the XY plane.
		self.setPlaneNormal(Position(-0.00000000001, 0.0, 1.0))
		#self.setPlaneNormal(Position(0.0, 0.0, 1.0))
		#self.setPlaneNormal(Position(-1, 0.0, 0.0001))
		#self.setPlaneNormal(Position(-0.287, 0.0, 0.9579))
		self.vtk_plane.SetTransform(self.transform)

	def setPlaneOrigin(self, position):
		"""
		Set the plane origin.
		@type position: L{Position <geo.Position>} object
		@param position: Plane origin
		"""

		self.vtk_plane.SetOrigin(position.getXCoor(), position.getYCoor(),
			position.getZCoor())

	def setPlaneNormal(self, position):
		"""
		Set the plance normal.
		@type position: L{Position <geo.Position>} object
		@param position: Plane normal 
		"""

		self.vtk_plane.SetNormal(position.getXCoor(), position.getYCoor(),
			position.getZCoor())

	def setCutter(self, component):
		"""
		Set up the cutter.
		@type component: String
		@param component: Component to be cut using the plane
		"""

		self.vtk_cutter.SetInput(component)
		self.vtk_cutter.SetCutFunction(self.vtk_plane)

	def setClipper(self, component):
		"""
		Set up the clipper.
		@type component: String
		@param component: Component to be clip using the plane
		"""

		self.vtk_clipper.SetInput(component)
		self.vtk_clipper.SetClipFunction(self.vtk_plane)
		#self.vtk_clipper.SetValue(0.0)
		# Use implicit function to clip and instead of input's scalar data
		#self.vtk_clipper.GenerateClipScalarsOn()
		# Generate the polygonal data that has been clipped away
		self.vtk_clipper.GenerateClippedOutputOn()
		# Specifies the clipping side of the plane.
		self.setInsideOutOn()

	def setValue(self, clipping_value):
		"""
		Set the clipping value.
		@type clipping_value: Number
		@param clipping_value: Clipping value of the implicit function or 
			scalar value
		"""

		self.vtk_clipper.SetValue(clipping_value)
	
	def setInsideOutOn(self):
		"""
		Set the clipping to inside out.
		"""

		self.vtk_clipper.InsideOutOn()

	def setInsideOutOff(self):
		"""
		Disable the inside out clipping.
		"""

		self.vtk_clipper.InsideOutOff()


#def Plane(object): ------------------------------------
"""
A plane in global coordinates
"""
pass

def Origin(Position):
    """
    The position of the origin
    """
    pass

def Direction(object):
    """
    A dirction in global coordinates
    """
    pass

def XAxis(Direction):
    """
    The direction of the x-axis
    """
    pass

def YAxis(Direction):
    """
    The direction of the y-axis
    """
    pass

def ZAxis(Direction):
    """
    The direction of the z-axis
    """
    pass


def YZPlane(Plane):
    """
    The YZ plane orthogonal to the x-axis
    """
    pass

def ZXPlane(Plane):
    """
    The ZX plane orthogonal to the y-axis
    """
    pass

def Sphere(object):
    """
    A sphere
    """
    pass

