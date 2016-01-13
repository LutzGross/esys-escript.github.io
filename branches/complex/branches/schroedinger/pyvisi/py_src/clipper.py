
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


import vtk

class Clipper:
	"""
	Class that defines a clipper.
	"""
	
	def __init__(self):
		"""
		Initialise the clipper.
		"""
		
		self.__vtk_clipper = vtk.vtkClipDataSet()

	def _setupClipper(self, object, plane):
		"""
		Setup the clipper.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the clipper
		@type plane: vtkPlane
		@param plane: Plane to clip the object 
		"""

		self.__object = object 
		# True only if a plane is used to perform clipping. False for 
		# scalar clipping.
		if(plane != None): 
			self.__plane = plane
		
		self.__setInput()
		# Due to this it's not possible to setInsideOutOff() when used with Rotation
		# self.setInsideOutOn()

	def __setInput(self):
		"""
		Set the input for the clipper.
		"""

		self.__vtk_clipper.SetInput(self.__object)

	def _setClipFunction(self):
		"""
		Set the clip function (using a plane) for the clipper.
		"""

		self.__vtk_clipper.SetClipFunction(self.__plane)

	def setInsideOutOn(self):
		"""
		Clip one side of the rendered object.		
		"""

		self.__vtk_clipper.InsideOutOn()

	def setInsideOutOff(self):
		"""
		Clips the other side of the rendered object.		
		"""

		self.__vtk_clipper.InsideOutOff()

	def setClipValue(self, value):
		"""
		Set the scalar clip value (intead of using a plane) for the clipper.

		@type value: Number
		@param value: Scalar clip value
		"""

		self.__vtk_clipper.SetValue(value)

	def _getClipperOutput(self):
		"""
		Return the output of the clipper.

		@rtype: vtkUnstructuredGrid
		@return: Unstructured grid
		"""

		return self.__vtk_clipper.GetOutput()

