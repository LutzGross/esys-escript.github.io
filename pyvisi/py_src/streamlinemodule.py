"""
@author: John NGUI
"""

import vtk

class StreamLineModule:
	"""
	Class that defines streamline.
	"""

	def __init__(self, object, source):
		"""
		Initialise the streamline.

		@type object: vtkUnstructuredGrid, etc 
		@param object: Input for the streamline
		@type source: vtkPolyData
		@param source: Source to generate the starting points
		"""

		self.__object = object
		self.__source = source
		self.__vtk_stream_line = vtk.vtkStreamLine()

		self.__setupStreamLineModule()

	def __setupStreamLineModule(self):
		"""
		Setup the streamline.
		"""

		self.__setInput()
		self.__setSource()
		# Default maximum propagation time is 100.
		self.setMaximumPropagationTime(100)
		# Default step length is 0.1
		self.setStepLength(0.1)
		# Default integration step length is 0.1
		self.setIntegrationStepLength(0.1)
		# Default integration is set to both directions.
		self.setIntegrationToBothDirections()
		# Default integrator is set to vtkRungeKutta4
		self.setIntegrator(vtk.vtkRungeKutta4())	
		self.__vtk_stream_line.Update()

	def __setInput(self):
		"""
		Set the input for the streamline.
		"""

		self.__vtk_stream_line.SetInput(self.__object)

	def __setSource(self):
		"""
		Set the source to generate the starting points for the streamlines.
		"""

		self.__vtk_stream_line.SetSource(self.__source)

	def setMaximumPropagationTime(self, time):
		"""
		Set the maximum length of the streamline expressed in elapsed time.

		@type time: Number
		@param time: Maximum length of the streamline expressed in elapsed time
		"""

		self.__vtk_stream_line.SetMaximumPropagationTime(time)

	def setStepLength(self, length):
		"""
		Set the length of the line segment expressed in elapsed time. A smaller
		value results in a smoother streamline (but is more expensive). Setting
		the step length usually goes hand-in-hand with setting the integration
		step length. Otherwise, errors such as "... can't compute normals" may
		arise. However, it does not usually apply the other way around.

		@type length: Number
		@param length: Length of the line segment expressed in elapsed time
		"""

		self.__vtk_stream_line.SetStepLength(length)	

	def setIntegrationStepLength(self, length):
		"""
		Set the integration step size expressed as a fraction of the size of 
		each cell. A smaller length results in a better image (but is more 
		expensive).

		@type length: Number
		@param length: Length of the integration step expressed as a fraction 
		"""

		self.__vtk_stream_line.SetIntegrationStepLength(length)

	def setIntegrationToBothDirections(self):
		"""
		Set the integration to occur both sides: forward (where the streamline
		goes) and backward (where the streamline came from).
		"""

		self.__vtk_stream_line.SetIntegrationDirectionToIntegrateBothDirections()

	def setIntegrator(self, integrator):
		"""
		Set the integrator to be used in the streamline calculation.

		@type integrator: vtkInitialValueProblemSolver
		@param integrator: Integrator type. i.e. vtkRungeKutta2, vtkRungeKutta4
		"""

		self.__vtk_stream_line.SetIntegrator(integrator)

	def _getOutput(self):
		"""
		Return the streamline.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__vtk_stream_line.GetOutput()
