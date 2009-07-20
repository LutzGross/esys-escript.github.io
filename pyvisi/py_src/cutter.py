
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

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class Cutter:
	"""
	Class that defines a cutter.
	"""

	def __init__(self):
		"""
		Initialise the cutter.
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Cutter is not running on more than one processor"
		self.__vtk_cutter = vtk.vtkCutter()

	def _setupCutter(self, object, plane):
		"""
		Setup the cutter.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the cutter
		@type plane: vtkPlane
		@param plane: Plane to cut the object
		"""

		self.__object = object
		self.__plane = plane

		self.__setInput()
		self.__setCutFunction()

	def __setInput(self):
		"""
		Set the input for the cutter.
		"""

		self.__vtk_cutter.SetInput(self.__object)

	def __setCutFunction(self):
		"""
		Set the cut function (using a plane).
		"""

		self.__vtk_cutter.SetCutFunction(self.__plane)

	def _getCutterOutput(self):
		"""
		Return the output of the cutter.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_cutter.GetOutput()
