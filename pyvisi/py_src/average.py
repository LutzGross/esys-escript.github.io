
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class CellDataToPointData:
	"""
	Class that defines a filter to convert cell data to point data by 
	averaging the data values of all cells using a particular point.
	"""

	def __init__(self, object):
		"""
		Initialise the cell to point data filter.

		:type object: vtkUnstructuredGrid, etc
		:param object: Input for the cell to point data filter
		"""
                if getMPISizeWorld()>1:
                    raise ValueError,"pyvisi.CellDataToPointData is not running on more than one processor."
		self.__object = object
		self.__vtk_cell_to_point = vtk.vtkCellDataToPointData()

		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the cell to point data filter.
		"""

		self.__vtk_cell_to_point.SetInput(self.__object)

	def _getCellToPointOutput(self):
		"""
		Return the output of the cell to point data filter.

		:rtype: vtkDataSet
		:return: Data set
		"""

		return self.__vtk_cell_to_point.GetOutput()

	
