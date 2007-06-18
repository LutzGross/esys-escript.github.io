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

class Sphere:
	"""
	Class that defines a sphere.
	"""

	def __init__(self):
		"""
		Initialise the sphere.
		"""

		self.__vtk_sphere = vtk.vtkSphereSource()

	def setThetaResolution(self, resolution):
		"""
		Set the theta resolution of the sphere.

		@type resolution: Number
		@param resolution: Theta resolution
		"""

		self.__vtk_sphere.SetThetaResolution(resolution)

	def setPhiResolution(self, resolution):
		"""
		Set the phi resolution of the sphere.

		@type resolution: Number
		@param resolution: Phi resolution
		"""

		self.__vtk_sphere.SetPhiResolution(resolution)

	def _getSphereOutput(self):
		"""
		Return the output of the sphere.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_sphere.GetOutput()

