"""
@author: John NGUI
"""

import vtk

class Probe:
	def __init__(self, object, source):
		self.__object = object
		self.__source = source
		self.__vtk_probe_filter = vtk.vtkProbeFilter()
		self.__setInput()
		self.__setSource()

	def __setInput(self):
		self.__vtk_probe_filter.SetInput(self.__source)

	def __setSource(self):
		self.__vtk_probe_filter.SetSource(self.__object)

	def _getOutput(self):
		return self.__vtk_probe_filter.GetOutput()

