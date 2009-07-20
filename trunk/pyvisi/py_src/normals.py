
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

class Normals:
	"""
	Class that defines normals. Normals are used to average the normals of 
	points in order to generate better sufaces (in the case of tensors, normals
	avoids the tensors from appearing black in color).
	"""

	def __init__(self):
		"""
		Initialise the normals.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.Normals is not running on more than one processor."
		self.__vtk_poly_data_normals = vtk.vtkPolyDataNormals()

	def _setupNormals(self, object):
		"""
		Setup the normals.

		@type object: vtkPolyData, etc
		@param object: Input for the normals
		"""

		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the normals.
		"""

		self.__vtk_poly_data_normals.SetInput(self.__object)

	def _getNormalsOutput(self):
		"""
		Return the output of the normals.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_poly_data_normals.GetOutput()
