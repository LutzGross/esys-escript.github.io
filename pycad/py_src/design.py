# $Id:$

"""
template for the Design which defines a regions and features
for a mesh generator.

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""


__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2007 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"


from primitives import Primitive, ReversePrimitive

class Design(object):
    """
    template for a design which defines the input for a mesh generator
    """
    def __init__(self,dim=3,element_size=1.,order=1,keep_files=False):
       """
       initializes a design 

       @param dim: patial dimension
       @element_size: global element size
       @order: element order
       @keep_files: flag to keep work files.
       """ 
       self.clearItems()
       self.setElementSize(element_size)
       self.setDim(dim)
       self.setElementOrder(order)
       if keep_files:
          self.setKeepFilesOn()
       else:
          self.setKeepFilesOff()
    def setDim(self,dim=3):
        """
        sets the spatial dimension
        """
        if not dim  in [1,2,3]: 
           raise ValueError("only dimension 1, 2, 3 are supported.")
        self.__dim=dim
    def getDim(self,dim=3):
        """
        returns the spatial dimension
        """
        return self.__dim
    def setElementOrder(self,order=1):
        """
        sets the element order
        """
        if not order in [1,2]:
           raise ValueError("only element orser 1 or 2 is supported.")
        self.__order=order
    def getElementOrder(self,order=1):
        """
        returns the element order
        """
        return self.__order
    def setElementSize(self,element_size=0.1):
        """
        set the global element size.
        """
        if element_size<=0.:
           raise ValueError("element size needs to be non--negative.")
        self.__element_size=element_size
    def getElementSize(self,element_size=1.):
        """
        returns the global element size.
        """
        return self.__element_size
    def setKeepFilesOn(self):
        """
        work files are kept at the end of the generation
        """
        self.__keep_files=True
    def setKeepFilesOff(self):
        """
        work files are deleted at the end of the generation
        """
        self.__keep_files=False
    def keepFiles(self):
        """
        returns True if work files are kept
        """
        return self.__keep_files
    def addItems(self,*items):
       """
       adds items to the design
       """
       for i in range(len(items)):
          if not isinstance(items[i],(Primitive, ReversePrimitive)):
             raise TypeError("%s-th argument is not a Primitive object"%i)
       for i in items:
          self.__items.append(i)
    def getItems(self):
        """
        returns a list of the items used in the design
        """
        return self.__items
    def clearItems(self):
        """
        resets the items in design
        """
        self.__items=[]
    def getAllPrimitives(self):
        """
        returns a list of all primitives used to create the design.
        each primitve appears once. The primitives are ordered by their
        order of generation
        """
        prims=set()
        for i in self.getItems(): prims|=set(i.getPrimitives())
        prims=list(prims)
        prims.sort()
        return prims

    def setOptions(self,**kwargs):
        """
        sets options of the mesh generator

        @note: this method is typically overwritten by a particular Design implementation
        """
        pass
    def getMeshHandler(self):
        """
        returns a handle to a mesh meshing the design

        @note: this method has to be overwritten by a particular Design implementation
        """
        raise NotImplementedError()

