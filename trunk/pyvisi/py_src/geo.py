"""
@author: John Ngui
@author: Lutz Gross
"""

class Position:
	"""
	Class that defines the positioning of components in the visualization.
	"""

	def __init__(self, x_coor, y_coor, z_coor):
		"""
		Initialize the x,y and z coordinates.	

		@type x_coor: Number
		@param x_coor: X coordinate in global position
		@type y_coor: Number
		@param y_coor: Y coordinate in global position
		@type z_coor: Number 
		@param z_coor: Z coordinate in global position
		"""
		
		self.x_coor = x_coor 
		self.y_coor = y_coor 
		self.z_coor = z_coor

	def getXCoor(self):
		"""	
		Return the x coordinate.

		@rtype: Number
		@return: X coordinate
		"""
		return self.x_coor

	def getYCoor(self):
		"""
		Return the y coordinate.

		@rtype: Number
		@return: Y coordiante
		"""

		return self.y_coor

	def getZCoor(self):
		"""
		Return the z coordinate

		@rtype: Number
		@return: Z coordinate
		"""

		return self.z_coor


from common import *

class Plane(Common):
	"""
	Class that performs cutting using a plane as its implicit function.
	"""

	def __init__(self, scene, data_collector, component):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@type component: String
		@param component: Component to be cut using the plane
		"""

		Common.__init__(self, scene, data_collector)
		self.vtk_plane = None
		self.vtk_cutter = None
		self.setPlane()
		self.setCutter(component)

		Common.setMapper(self, "self.vtk_cutter.GetOutput()")
		Common.setActor(self)
		Common.addActor(self)

	def setPlane(self):
		"""
		Setup the plane.
		"""

		self.vtk_plane = vtk.vtkPlane()
		# Default origin
		#self.vtk_plane.SetOrigin(
		#self.data_collector.getReader().GetOutput().GetCenter())
		self.vtk_plane.SetOrigin(0,0,0)
		# Default normal
		self.vtk_plane.SetNormal(-1.2, 0.0, 0.9)


	def setPlaneOrigin(self, x_coor, y_coor, z_coor):
		"""
		Set the plane origin.

		@type x_coor: Number
		@param x_coor: X coordinate in global position
		@type y_coor: Number
		@param y_coor: Y coordinate in global position
		@type z_coor: Number
		@param z_coor: Z coordinate in global position
		"""

		self.vtk_plane.SetOrigin(x_coor, y_coor, z_coor)

	def setPlaneNormal(self, x_coor, y_coor, z_coor):
		"""
		Set the plance normal.

		@type x_coor: Number
		@param x_coor: X coordinate in global position
		@type y_coor: Number
		@param y_coor: Y coordinate in global position
		@type z_coor: Number
		@param z_coor: Z coordinate in global position
		"""

		self.vtk_plane.SetNormal(x_coor, y_coor, z_coor)

	def setCutter(self, component):
		"""
		Setup the cutter

		@type component: String
		@param component: Component to be cut using the plane
		"""

		self.vtk_cutter = vtk.vtkCutter()
		eval("self.vtk_cutter.SetInput(%s)" % component)
		self.vtk_cutter.SetCutFunction(self.vtk_plane)

#def Plane(object):
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


def XYPlane(Plane):
    """
    The XY plane orthogonal to the z-axis
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

