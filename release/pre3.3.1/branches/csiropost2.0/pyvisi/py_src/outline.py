
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

class Outline:
	"""
	Class that defines an outline.
	"""

	def __init__(self, object):
		"""
		Initialise the outline.

		@type object: vtkUnstructuredGrid, etc
		@param object: Data source to the outline
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.Outline is not running on more than one processor."
		self.__object = object
		self.__vtk_outline = vtk.vtkOutlineFilter()
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the outline.
		"""

		self.__vtk_outline.SetInput(self.__object)
	
	def _getOutlineOutput(self):
		"""
		Return the output of the outline.

		@rtype: vtkPolyData
		@return: Polyognal data
		"""

		return self.__vtk_outline.GetOutput()
	
