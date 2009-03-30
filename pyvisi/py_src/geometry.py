
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

class Geometry:
	"""
	Class that extracts geometry from data and convert it to polygonal type.	
	"""

	def __init__(self, object):	
		"""
		Initialise the geometry filter.
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Geometry is not running on more than one processor"
		self.__vtk_geometry_filter = vtk.vtkGeometryFilter()
		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the geometry filter
		"""

		self.__vtk_geometry_filter.SetInput(self.__object)

	def _getGeometryOutput(self):
		"""
		Return the output of the rotation.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""
		return self.__vtk_geometry_filter.GetOutput()
