"""

classes dealing with data for the visualization

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane, L. Gross"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision:$"
__date__="$Date:$"

import vtk

class DataCollector:

	def __init__(self, outline = True):
		self.outline = True
		self.file_name = None
		self.vtk_xml_reader = None

	def setSource(self, file_name):
		self.file_name = file_name
		self.vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
		self.vtk_xml_reader.SetFileName(self.file_name)

	def getReader(self):
		return self.vtk_xml_reader


