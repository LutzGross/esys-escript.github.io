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
	
