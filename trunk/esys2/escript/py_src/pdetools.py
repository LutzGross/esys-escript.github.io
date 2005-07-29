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
    self.__pde.setLumping(fast)
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


class Location:
     """Location provides a factory to access the values of data objects at a given spatial coordinate x.
        In fact, a Location object finds the sample in the set of samples of a given function space or domain which is closest
        to the given point x. """

     def __init__(self,x,where):
       """initializes a Location to access values in Data objects on the Doamin or FunctionSpace where for the sample point which
          closest to the given point x"""
       if isinstance(where,FunctionSpace):
          self.__where=where
       else:
          self.__where=ContinuousFunction(where)
       self.__id=length(x-self.__where.getX()).minarg()

     def getValue(self,data):
        """returns the value of data at the Location at a numarray object"""
        if isinstance(data,Data):
           return data
        else:
           if not data.getFunctionSpace()==self.getFunctionSpace():
             raise ValueError,"function space of data obejct does not match function space of Location"
           else:
             return data.getValue(self.getId())
     def getX(self):
        """returns the exact coordinates of the Location"""
        return self.getValue(self.getFunctionSpace().getX())

     def getId(self):
        """returns the identifier of the location"""
        return self.__id

     def getFunctionSpace(self):
        """returns the function space of the Location"""
        return self.__function_space

     def __str__(self):
       """returns the coordinates of the Location as a string"""
       return str(self.getX())

def testProjector(domain):
      """runs a few test of the Projector factory and returns the largest error plus a description of the test this error occured"""
      error_max=0.
      error_text=""
      x=ContinuousFunction(domain).getX()
      for f in [True,False]:
         p=Projector(domain,reduce=False,fast=f)
         for r in range(5):
            text="range %s , fast=%s"%(r,f)
            if r==0:
               td_ref=x[0]
            elif r==1:
               td_ref=x
            elif r==2:
               td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
            elif r==3:
               td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
            elif r==3:
               td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
            td=p(td_ref.interpolate(Function(domain)))
            er=Lsup(td-td_ref)/Lsup(td_ref)
            print text," error = ",er
            if er>error_max:
                 error_max=er
                 error_text=text
      return error_max,error_text


# this should be removed later
if __name__=="__main__":
    from esys.finley import Rectangle
    txt=testProjector(Rectangle(56,61))
    print "test Projector: ",txt



