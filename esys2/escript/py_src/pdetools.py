# $Id$

"""
provides a some tools related to PDEs currently includes:

     Projector - to project a discontinuous


"""

from escript import *
from linearPDEs import LinearPDE
import numarray

class Projector:
  """The Projector is a factory which projects a discontiuous function onto a continuous function on the a given domain"""
  def __init__(self, domain, reduce = True, fast=True):
    """
    @brief Create a continuous function space projector for a domain.

    @param domain Domain of the projection.
    @param reduce Flag to reduce projection order (default is True)
    @param fast Flag to use a fast method based on matrix lumping (default is true)
    """
    self.__pde = LinearPDE(domain)
    if fast:
      self.__pde.setSolverMethod(LinearPDE.LUMPING)
    self.__pde.setSymmetryOn()
    self.__pde.setReducedOrderTo(reduce)
    self.__pde.setValue(D = 1.)
    return

  def __del__(self):
    return

  def __call__(self, input_data):
    """
    @brief projects input_data onto a continuous function

    @param input_data  The input_data to be projected.
    """
    out=Data(0.,input_data.getShape(),what=ContinuousFunction(self.__pde.getDomain()))
    if input_data.getRank()==0:
        self.__pde.setValue(Y = input_data)
        out=self.__pde.getSolution()
    elif input_data.getRank()==1:
        for i0 in range(input_data.getShape()[0]):
           self.__pde.setValue(Y = input_data[i0])
           out[i0]=self.__pde.getSolution()
    elif input_data.getRank()==2:
        for i0 in range(input_data.getShape()[0]):
           for i1 in range(input_data.getShape()[1]):
              self.__pde.setValue(Y = input_data[i0,i1])
              out[i0,i1]=self.__pde.getSolution()
    elif input_data.getRank()==3:
        for i0 in range(input_data.getShape()[0]):
           for i1 in range(input_data.getShape()[1]):
              for i2 in range(input_data.getShape()[2]):
                 self.__pde.setValue(Y = input_data[i0,i1,i2])
                 out[i0,i1,i2]=self.__pde.getSolution()
    else:
        for i0 in range(input_data.getShape()[0]):
           for i1 in range(input_data.getShape()[1]):
              for i2 in range(input_data.getShape()[2]):
                 for i3 in range(input_data.getShape()[3]):
                    self.__pde.setValue(Y = input_data[i0,i1,i2,i3])
                    out[i0,i1,i2,i3]=self.__pde.getSolution()
    return out


class Locator:
     """

        Locator provides a  access the values of data objects at a given spatial coordinate x.
        In fact, a Locator object finds the sample in the set of samples of a given function space or domain where which is closest
        to the given point x. 

     """

     def __init__(self,where,x=numarray.zeros((3,))):
       """initializes a Locator to access values in Data objects on the Doamin or FunctionSpace where for the sample point which
          closest to the given point x"""
       if isinstance(where,FunctionSpace):
          self.__function_space=where
       else:
          self.__function_space=ContinuousFunction(where)
       self.__id=length(x[:self.__function_space.getDim()]-self.__function_space.getX()).mindp()

     def __str__(self):
       """returns the coordinates of the Locator as a string"""
       return "<Locator %s>"%str(self.getX())

     def getFunctionSpace(self):
        """returns the function space of the Locator"""
        return self.__function_space

     def getId(self):
        """returns the identifier of the location"""
        return self.__id

     def getX(self):
        """returns the exact coordinates of the Locator"""
        return self(self.getFunctionSpace().getX())

     def __call__(self,data):
        """returns the value of data at the Locator of a Data object otherwise the object is returned."""
        return self.getValue(data)

     def getValue(self,data):
        """returns the value of data at the Locator if data is a Data object otherwise the object is returned."""
        if isinstance(data,Data):
           if data.getFunctionSpace()==self.getFunctionSpace():
             out=data.convertToNumArrayFromDPNo(self.getId()[0],self.getId()[1])
           else:
             out=data.interpolate(self.getFunctionSpace()).convertToNumArrayFromDPNo(self.getId()[0],self.getId()[1])
           if data.getRank()==0:
              return out[0]
           else:
              return out
        else:
           return data
