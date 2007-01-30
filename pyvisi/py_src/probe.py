"""
@author: John NGUI
"""

import vtk

class Probe:
	"""
	Class that defines the probe. The class sample data values at 
	specified point locations.
	"""

	def __init__(self, object, source):
		"""
		Initialise the probe.

		@type object: vtkUnstructuredGrid, etc
		@param object: Inpur for the probe 
		@type source: vtkDataSet (i.e. vtkStructuredPoints)
		@param source: Source for the probe
		"""

		self.__object = object
		self.__source = source
		self.__vtk_probe_filter = vtk.vtkProbeFilter()
		self.__setInput()
		self.__setSource()

	def __setInput(self):
		"""
		Set the input for the probe.
		"""

		self.__vtk_probe_filter.SetInput(self.__source)

	def __setSource(self):
		"""
		Set the source for the probe.
		"""

		self.__vtk_probe_filter.SetSource(self.__object)

	def _getOutput(self):
		"""
		Return the probe.

		@rtype: vtkDataSet
		@return: Data set
		"""

		return self.__vtk_probe_filter.GetOutput()

