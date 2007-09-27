#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
@author: John NGUI
"""

import vtk

class Probe:
	"""
	Class that defines the probe. This class samples data values at 
	the specified point locations.
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

		self.__setupProbe()

	def __setupProbe(self):
		"""
		Setup the probe.
		"""

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
		Return the output for the probe.

		@rtype: vtkDataSet
		@return: Data set
		"""

		return self.__vtk_probe_filter.GetOutput()

