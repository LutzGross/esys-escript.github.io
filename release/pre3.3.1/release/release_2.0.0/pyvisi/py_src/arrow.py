
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

class Arrow2D:
	"""
	Class that defines 2D arrows.
	"""
	
	def __init__(self):
		"""
		Initialise the 2D arrows.
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Arrow2D is not running on more than one processor"
		self.__vtk_arrow2D = vtk.vtkGlyphSource2D()
		self.__setupArrow2D()

	def __setupArrow2D(self):
		"""
		Setup the 2D arrows.
		"""

		# Use arrows instead of cone or sphere.
		self.__vtk_arrow2D.SetGlyphTypeToArrow()
		# Fill the inside of the arrows.
		self.__vtk_arrow2D.SetFilled(0)
	
	def _getArrow2DOutput(self):
		"""
		Return the output of the 2D arrows.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""
	
		return self.__vtk_arrow2D.GetOutput()


###############################################################################


class Arrow3D:
	"""
	Class that defines 3D arrows.
	"""

	def __init__(self):	
		"""
		Initialise the 3D arrows.
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.Arrow3D is not running on more than one processor"
		self.__vtk_arrow3D = vtk.vtkArrowSource()
		
	def _getArrow3DOutput(self):
		"""
		Return the output of the 3D arrows.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__vtk_arrow3D.GetOutput()

