
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


from constant import WarpMode
from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class Warp:
	"""
	Class that defines the deformation of a scalar field.
	"""

	def __init__(self,  warp_mode):
		"""
		Initialise the warp scalar/vector.

		@type warp_mode: L{WarpMode <constant.WarpMode>} constant
		@param warp_mode: Mode in which to deform the data
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Warp is not running on more than one processor."
		if(warp_mode == WarpMode.SCALAR): # Deform data with scalar data.
			self.__vtk_warp = vtk.vtkWarpScalar()
		elif(warp_mode == WarpMode.VECTOR): # Deform data with vector data.
			self.__vtk_warp = vtk.vtkWarpVector()

	def _setupWarp(self, object):
		"""
		Setup the warp.

		@type object: vtkPolyData, etc
		@param object: Input for the warp scalar or warp vector
		"""

		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the warp scalar/vector.
		"""

		self.__vtk_warp.SetInput(self.__object)

	def setScaleFactor(self, scale_factor):
		"""
		Set the displacement scale factor.

		@type scale_factor: Number
		@param scale_factor: Scale factor for the displacement
		"""

		self.__vtk_warp.SetScaleFactor(scale_factor)

	def _getWarpOutput(self):
		"""
		Return the output of the deformed data. 
		
		@rtype: vtkPointSet
		@return: PointSet data
		"""
		
		return self.__vtk_warp.GetOutput()

