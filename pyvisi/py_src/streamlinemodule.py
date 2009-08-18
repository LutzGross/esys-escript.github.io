
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class StreamLineModule:
	"""
	Class that defines the streamline module.
	"""

	def __init__(self):
		"""
		Initialise the streamline module.
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.StreamLineModule is not running on more than one processor."
		self.__vtk_stream_line = vtk.vtkStreamLine()

	def _setupStreamLineModule(self, object, source):
		"""
		Setup the streamline.

		:type object: vtkUnstructuredGrid, etc 
		:param object: Input for the streamline
		:type source: vtkPolyData
		:param source: Source to generate the starting points
		"""

		self.__object = object
		self.__source = source

		self.__setInput()
		self.__setSource()
		# Default maximum propagation time is 100.
		self.setMaximumPropagationTime(100)
		# Default step length is 0.01
		self.setStepLength(0.01)
		# Default integration step length is 0.01
		self.setIntegrationStepLength(0.01)
		# Default integration is set to both directions.
		self.setIntegrationToBothDirections()
		# Default integrator is set to vtkRungeKutta4
		self.setIntegrator(vtk.vtkRungeKutta4())	

	def __setInput(self):
		"""
		Set the input for the streamline.
		"""

		self.__vtk_stream_line.SetInput(self.__object)

	def __setSource(self):
		"""
		Set the source to generate the starting points for the streamline.
		"""

		self.__vtk_stream_line.SetSource(self.__source)

	def setMaximumPropagationTime(self, time):
		"""
		Set the maximum length of the streamline expressed in elapsed time.

		:type time: Number
		:param time: Maximum length of the streamline expressed in elapsed time
		"""

		self.__vtk_stream_line.SetMaximumPropagationTime(time)

	def setStepLength(self, length):
		"""
		Set the length of the streamline segment expressed in elapsed time. 
		A smaller value results in a smoother streamline 
		(but is more expensive). Setting the step length usually goes 
		hand-in-hand with setting the integration step length. Otherwise, 
		errors such as "... can't compute normals" may arise. If such an 
		error occurs try changing the values. 

		:type length: Number
		:param length: Length of the streamline segment expressed in 
				elapsed time
		"""

		self.__vtk_stream_line.SetStepLength(length)	

	def setIntegrationStepLength(self, length):
		"""
		Set the integration step size expressed as a fraction of the size of 
		each cell. A smaller length results in a better image (but is more 
		expensive).

		:type length: Number
		:param length: Length of the integration step expressed as a fraction 
				of the size of each cell
		"""

		self.__vtk_stream_line.SetIntegrationStepLength(length)

	def setIntegrationToBothDirections(self):
		"""
		Set the integration to occur both sides: forward (where the streamline
		goes) and backward (where the streamline came from).
		"""

		self.__vtk_stream_line.\
				SetIntegrationDirectionToIntegrateBothDirections()

	def setIntegrator(self, integrator):
		"""
		Set the integrator to be used in the streamline calculation. 

		:type integrator: vtkInitialValueProblemSolver
		:param integrator: Integrator type. i.e. vtkRungeKutta2, vtkRungeKutta4
		"""

		self.__vtk_stream_line.SetIntegrator(integrator)
	
	def _setSpeedScalarsOn(self):
		"""
		Turn on the creation of scalar data from velocity magnitude. 
		"""

		self.__vtk_stream_line.SpeedScalarsOn()

	def _setSpeedScalarsOff(self):
		"""
		Turn off the creation of scalar data from velocity magnitude. If
		input dataset has scalars, input dataset scalars are used.
		"""

		self.__vtk_stream_line.SpeedScalarsOff()

	def _getStreamLineModuleOutput(self):
		"""
		Return the output of the streamline.

		:rtype: vtkPolyData
		:return Polygonal data
		"""

		return self.__vtk_stream_line.GetOutput()
