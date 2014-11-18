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
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


import vtk

class Rotation:
	"""
	Class that sweeps 2D data around the z-axis to create a 3D looking effect.
	"""

	def __init__(self):
		"""
		Initialise the rotation.
		"""

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
		