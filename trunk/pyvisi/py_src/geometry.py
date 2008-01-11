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

class Geometry:
	"""
	Class that extracts geometry from data and convert it to polygonal type.	
	"""

	def __init__(self, object):	
		"""
		Initialise the geometry filter.
		"""

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
