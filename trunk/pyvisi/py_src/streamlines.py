"""
@author: John Ngui
@author: Lutz Gross
"""

import vtk
from common import Common
from geo import Position

class StreamLines(Common):
	"""
	Class that shows the path of particles (within a specified cloud of points)
	in a vector field.
	"""

	def __init__(self, scene, data_collector, lut = None):
		"""
		@type scene: L{Scene <scene.Scene>} object
		@param scene: Scene in which components are to be added to
		@type data_collector: L{DataCollector <datacollector.DataCollector>}
			object
		@param data_collector: Source of data for visualization
		@type lut: L{BlueToRed <colormap.BlueToRed>} or
			L{RedToBlue <colormap.RedToBlue>} object
		@param lut: Lookup table to be used by the mapper
		"""

		Common.__init__(self, scene, data_collector)
		# Create a spherical cloud of points.
		self.vtk_point_source = vtk.vtkPointSource()
		# Create streamlines
		self.vtk_stream_lines = vtk.vtkStreamLine()
		# Create tube to be wrapped around the streamlines.
		self.vtk_tube = vtk.vtkTubeFilter()

		self.setSource()
		self.setStreamLines()
		self.setTube()
		
		Common.setMapperInput(self, self.vtk_tube.GetOutput(), lut)
		Common.setActorInput(self)
		Common.addActor(self)

	def setSource(self):
		"""
		Set the default parameters for the spherical cloud of points.
		"""

		self.setCloudRadius(0.2)
		self.setCenter(Position(0.2, 2.1, 0.5))	
		self.setNumberOfPoints(8)

	def setStreamLines(self):
		"""
		Set the default parameters for the streamlines.
		"""

		self.vtk_stream_lines.SetInput(
			self.data_collector.getReader().GetOutput())
		self.vtk_stream_lines.SetSource(self.vtk_point_source.GetOutput())
		self.setMaximumPropagationTime(500)
		self.setLineSize(0.2)
		self.setAccuracy(0.05)
		
		self.setIntegrationToBothDirections()
		self.vtk_stream_lines.SetIntegrator(vtk.vtkRungeKutta4())


	def setTube(self):
		"""
		Set the default parameters for the tube.
		"""

		self.vtk_tube.SetInput(self.vtk_stream_lines.GetOutput())
		self.setTubeRadius(0.01)
		self.setNumberOfSides(12)
		self.setVaryRadiusByVector()
		
	def setCloudRadius(self, radius):
		"""
		Set the radius for the cloud of points.
		@type radius: Number
		@param radius: Radius for the cloud of points
		"""

		self.vtk_point_source.SetRadius(radius)

	def setCenter(self, position):
		"""
		Set the center for the cloud of points.
		@type position: L{Position <geo.Position>} object
		@param position: Center for the cloud of points
		"""

		self.vtk_point_source.SetCenter(position.getXCoor(), 
			position.getYCoor(), position.getZCoor())

	def setNumberOfPoints(self, points):
		"""
		Set the number of points to generate for the cloud of points.
		@type points: Number
		@param points: Number of points to generate for the cloud of points
		"""

		self.vtk_point_source.SetNumberOfPoints(points)

	def setMaximumPropagationTime(self, time):
		"""
		Set the maximum length for the streamlines in units of time.
		@type time: Number
		@param time: Units of time specifying the maximum length of streamlines 
		"""

		self.vtk_stream_lines.SetMaximumPropagationTime(time)

	def setStreamLinesSize(self, stream_lines_size):
		"""
		Set the size of the streamlines.
		@type stream_lines_size: Number
		@param stream_lines_size: Size of the streamlines
		"""

		self.vtk_stream_lines.SetStepLength(stream_lines_size)

	def setAccuracy(self, accuracy):
		"""
		Set the level of accuracy for the streamlines. The smaller the value 
		the more accurate it is.
		@type accuracy: Number
		@param accuracy: Accuracy for the streamlines (between 0 to 1)
		"""

		self.vtk_stream_lines.SetIntegrationStepLength(accuracy)

	def setIntegrationToBothDirections(self):
		"""
		Set the integration to occur in both directions.
		"""		

		self.vtk_stream_lines.SetIntegrationDirectionToIntegrateBothDirections()

	def setTubeRadius(self, radius):
		"""
		Set the minimum radius of the tube.
		@type radius: Number
		@param radius: Minimum radius of the tube
		"""

		self.vtk_tube.SetRadius(radius)

	def setNumberOfSides(self, sides):
		"""
		Set the number of sides for the tube.
		@type sides: Number
		@param sides: Number of sides for the tube (mimimum is 3)
		"""

		self.vtk_tube.SetNumberOfSides(sides)

	def setVaryRadiusByVector(self):
		"""
		Set the variation of the tube radius with vector data.
		"""
		
		self.vtk_tube.SetVaryRadiusToVaryRadiusByVector()
