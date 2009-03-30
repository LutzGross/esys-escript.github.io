
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


import vtk
from esys.escript import getMPISizeWorld

class Rotation:
	"""
	Class that sweeps 2D data around the z-axis to create a 3D looking effect.
	"""

	def __init__(self):
		"""
		Initialise the rotation.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.Rotation is not running on more than one processor."
		self.__vtk_rotational_extrusion_filter = \
				vtk.vtkRotationalExtrusionFilter()

	def _setupRotationExtrusionFilter(self, object):
		"""
		Setup the rotation.

		@type object: vtkUnstructuredGrid, etc
		@param object: Input for the clipper
		"""

		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the rotation.
		"""

		self.__vtk_rotational_extrusion_filter.SetInput(self.__object)
		self.__vtk_rotational_extrusion_filter.Update()
		
	def setResolution(self, resolution):
		"""
		Set the resolution of the sweep for the rotation, which controls the 
		number of intermediate points

		@type resolution: Number
		@param resolution: Resolution of the rotation
		"""

		self.__vtk_rotational_extrusion_filter.SetResolution(resolution)

	def setAngle(self, angle):
		"""
		Set the angle of rotation.

		@type angle: Number
		@param angle: Angle of rotation
		"""

		self.__vtk_rotational_extrusion_filter.SetAngle(angle)
	
	def _getRotationOutput(self):
		"""
		Return the output of the rotation.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_rotational_extrusion_filter.GetOutput()
		
