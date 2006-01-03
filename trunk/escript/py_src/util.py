# $Id$
#
#      COPYRIGHT ACcESS 2004 -  All Rights Reserved
#
#   This software is the property of ACcESS.  No part of this code
#   may be copied in any form or by any means without the expressed written
#   consent of ACcESS.  Copying, use or modification of this software
#   by any unauthorised person is illegal unless that
#   person has a software license agreement with ACcESS.
#

"""
Utility functions for escript

@remark:  This module is under construction and is still tested!!!

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""
                                                                                                                                                                                                     
__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"


import math
import numarray
import escript
import os

# missing tests:

# def pokeShape(arg):
# def pokeDim(arg):
# def commonShape(arg0,arg1):
# def commonDim(*args):
# def testForZero(arg):
# def matchType(arg0=0.,arg1=0.):
# def matchShape(arg0,arg1):

# def transpose(arg,axis=None):
# def trace(arg,axis0=0,axis1=1):
# def reorderComponents(arg,index):

# def integrate(arg,where=None):
# def interpolate(arg,where):
# def div(arg,where=None):
# def grad(arg,where=None):

#
# slicing: get
#          set
#
# and derivatives

#=========================================================
#   some helpers:
#=========================================================
def saveVTK(filename,domain=None,**data):
    """
    writes a L{Data} objects into a files using the the VTK XML file format.

    Example:

       tmp=Scalar(..)
       v=Vector(..)
       saveVTK("solution.xml",temperature=tmp,velovity=v)

    tmp and v are written into "solution.xml" where tmp is named "temperature" and v is named "velovity"

    @param filename: file name of the output file
    @type filename: C{str}
    @param domain: domain of the L{Data} object. If not specified, the domain of the given L{Data} objects is used.
    @type domain: L{escript.Domain}
    @keyword <name>: writes the assigned value to the VTK file using <name> as identifier.
    @type <name>: L{Data} object.
    @note: The data objects have to be defined on the same domain. They may not be in the same L{FunctionSpace} but one cannot expect that all L{FunctionSpace} can be mixed. Typically, data on the boundary and data on the interior cannot be mixed.
    """
    if domain==None:
       for i in data.keys():
          if not data[i].isEmpty(): domain=data[i].getFunctionSpace().getDomain()
    if domain==None:
        raise ValueError,"no domain detected."
    else:
        domain.saveVTK(filename,data)

def saveDX(filename,domain=None,**data):
    """
    writes a L{Data} objects into a files using the the DX file format.

    Example:

       tmp=Scalar(..)
       v=Vector(..)
       saveDX("solution.dx",temperature=tmp,velovity=v)

    tmp and v are written into "solution.dx" where tmp is named "temperature" and v is named "velovity".

    @param filename: file name of the output file
    @type filename: C{str}
    @param domain: domain of the L{Data} object. If not specified, the domain of the given L{Data} objects is used.
    @type domain: L{escript.Domain}
    @keyword <name>: writes the assigned value to the DX file using <name> as identifier. The identifier can be used to select the data set when data are imported into DX.
    @type <name>: L{Data} object.
    @note: The data objects have to be defined on the same domain. They may not be in the same L{FunctionSpace} but one cannot expect that all L{FunctionSpace} can be mixed. Typically, data on the boundary and data on the interior cannot be mixed.
    """
    if domain==None:
       for i in data.keys():
          if not data[i].isEmpty(): domain=data[i].getFunctionSpace().getDomain()
    if domain==None:
        raise ValueError,"no domain detected."
    else:
        domain.saveDX(filename,data)

def kronecker(d=3):
   """
   return the kronecker S{delta}-symbol

   @param d: dimension or an object that has the C{getDim} method defining the dimension
   @type d: C{int} or any object with a C{getDim} method
   @return: the object u of rank 2 with M{u[i,j]=1} for M{i=j} and M{u[i,j]=0} otherwise
   @rtype d: L{numarray.NumArray} of rank 2.
   @remark: the function is identical L{identity}
   """
   return identityTensor(d)

def identity(shape=()):
   """
   return the shape x shape identity tensor

   @param shape: input shape for the identity tensor
   @type shape: C{tuple} of C{int}
   @return: array of shape shape x shape with M{u[i,k]=1} for M{i=k} and M{u[i,k]=0} otherwise for len(shape)=1 and, for len(shape)=2: M{u[i,j,k,l]=1} for M{i=k and j=l} and M{u[i,j,k,l]=0} otherwise.
   @rtype: L{numarray.NumArray} of rank 1, rankk 2 or rank 4.
   @raise ValueError: if len(shape)>2.
   """
   if len(shape)>0:
      out=numarray.zeros(shape+shape,numarray.Float)
      if len(shape)==1:
          for i0 in range(shape[0]):
             out[i0,i0]=1.

      elif len(shape)==2:
          for i0 in range(shape[0]):
             for i1 in range(shape[1]):
                out[i0,i1,i0,i1]=1.
      else:
          raise ValueError,"identity: length of shape is restricted to 2."
   else:
      out=1.
   return out

def identityTensor(d=3):
   """
   return the dxd identity matrix

   @param d: dimension or an object that has the C{getDim} method defining the dimension
   @type d: C{int} or any object with a C{getDim} method
   @return: the object u of rank 2 with M{u[i,j]=1} for M{i=j} and M{u[i,j]=0} otherwise
   @rtype: L{numarray.NumArray} of rank 2.
   """
   if hasattr(d,"getDim"):
      d=d.getDim()
   return identity(shape=(d,))

def identityTensor4(d=3):
   """
   return the dxdxdxd identity tensor

   @param d: dimension or an object that has the C{getDim} method defining the dimension
   @type d: C{int} or any object with a C{getDim} method
   @return: the object u of rank 4 with M{u[i,j,k,l]=1} for M{i=k and j=l} and M{u[i,j,k,l]=0} otherwise
   @rtype: L{numarray.NumArray} of rank 4.
   """
   if hasattr(d,"getDim"):
      d=d.getDim()
   return identity((d,d))

def unitVector(i=0,d=3):
   """
   return a unit vector u of dimension d with nonzero index i:

   @param i: index 
   @type i: C{int} 
   @param d: dimension or an object that has the C{getDim} method defining the dimension
   @type d: C{int} or any object with a C{getDim} method
   @return: the object u of rank 1 with M{u[j]=1} for M{j=i} and M{u[i]=0} otherwise
   @rtype: L{numarray.NumArray} of rank 1.
   """
   return kronecker(d)[i]

#=========================================================================
#   global reduction operations (these functions have no symbolic version)
#=========================================================================
def Lsup(arg):
    """
    returns the Lsup-norm of argument arg. This is the maximum absolute value over all data points. 
    This function is equivalent to sup(abs(arg)).

    @param arg: argument
    @type arg: C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
    @return: maximum value of the absolute value of arg over all components and all data points
    @rtype: C{float}
    @raise TypeError: if type of arg cannot be processed
    """
    if isinstance(arg,numarray.NumArray):
        return sup(abs(arg))
    elif isinstance(arg,escript.Data):
        return arg._Lsup()
    elif isinstance(arg,float):
        return abs(arg)
    elif isinstance(arg,int):
        return abs(float(arg))
    else:
      raise TypeError,"Lsup: Unknown argument type."

def sup(arg):
    """
    returns the maximum value over all data points. 

    @param arg: argument
    @type arg: C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
    @return: maximum value of arg over all components and all data points
    @rtype: C{float}
    @raise TypeError: if type of arg cannot be processed
    """
    if isinstance(arg,numarray.NumArray):
        return arg.max()
    elif isinstance(arg,escript.Data):
        return arg._sup()
    elif isinstance(arg,float):
        return arg
    elif isinstance(arg,int):
        return float(arg)
    else:
      raise TypeError,"sup: Unknown argument type."

def inf(arg):
    """
    returns the maximum value over all data points. 

    @param arg: argument
    @type arg: C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
    @return : minimum value of arg over all components and all data points
    @rtype: C{float}
    @raise TypeError: if type of arg cannot be processed
    """
    if isinstance(arg,numarray.NumArray):
        return arg.min()
    elif isinstance(arg,escript.Data):
        return arg._inf()
    elif isinstance(arg,float):
        return arg
    elif isinstance(arg,int):
        return float(arg)
    else:
      raise TypeError,"inf: Unknown argument type."


#=========================================================================
#   some little helpers
#=========================================================================
def pokeShape(arg):
    """
    identifies the shape of its argument

    @param arg: a given object 
    @type arg: L{numarray.NumArray},L{escript.Data},C{float}, C{int}, C{Symbol}
    @return: the shape of the argument
    @rtype: C{tuple} of C{int}
    @raise TypeError: if type of arg cannot be processed
    """

    if isinstance(arg,numarray.NumArray):
        return arg.shape
    elif isinstance(arg,escript.Data):
        return arg.getShape()
    elif isinstance(arg,float):
        return ()
    elif isinstance(arg,int):
        return ()
    elif isinstance(arg,Symbol):
        return arg.getShape()
    else:
      raise TypeError,"pokeShape: cannot identify shape"

def pokeDim(arg):
    """
    identifies the spatial dimension of its argument

    @param arg: a given object 
    @type arg: any
    @return: the spatial dimension of the argument, if available, or C{None}
    @rtype: C{int} or C{None}
    """

    if isinstance(arg,escript.Data):
        return arg.getFunctionSpace().getDim()
    elif isinstance(arg,Symbol):
        return arg.getDim()
    else:
        return None

def commonShape(arg0,arg1):
    """
    returns a shape to which arg0 can be extendent from the right and arg1 can be extended from the left.

    @param arg0: an object with a shape (see L{pokeShape})
    @param arg1: an object with a shape (see L{pokeShape})
    @return: the shape of arg0 or arg1 such that the left port equals the shape of arg0 and the right end equals the shape of arg1.
    @rtype: C{tuple} of C{int}
    @raise ValueError: if no shape can be found.
    """
    sh0=pokeShape(arg0)
    sh1=pokeShape(arg1)
    if len(sh0)<len(sh1):
       if not sh0==sh1[:len(sh0)]:
             raise ValueError,"argument 0 cannot be extended to the shape of argument 1"
       return sh1
    elif len(sh0)>len(sh1):
       if not sh1==sh0[:len(sh1)]:
             raise ValueError,"argument 1 cannot be extended to the shape of argument 0"
       return sh0
    else:
       if not sh0==sh1:
             raise ValueError,"argument 1 and argument 0 have not the same shape."
       return sh0

def commonDim(*args):
    """
    identifies, if possible, the spatial dimension across a set of objects which may or my not have a spatial dimension.

    @param *args: given objects
    @return: the spatial dimension of the objects with identifiable dimension (see L{pokeDim}). If none the objects has 
             a spatial dimension C{None} is returned.
    @rtype: C{int} or C{None}
    @raise ValueError: if the objects with identifiable dimension don't have the same spatial dimension.
    """
    out=None
    for a in args:
       d=pokeDim(a)
       if not out==None:
          if not (d==None or out==d):
             raise ValueError,"dimension of arguments don't match"
       else:
          out=d
    return out

def testForZero(arg):
    """
    test the argument for being identical to Zero

    @param arg: a given object 
    @type arg: typically L{numarray.NumArray},L{escript.Data},C{float}, C{int}
    @return : True if the argument is identical to zero.
    @rtype : C{bool}
    """
    if isinstance(arg,numarray.NumArray):
       return not Lsup(arg)>0.
    elif isinstance(arg,escript.Data):
       return False
    elif isinstance(arg,float):
       return not Lsup(arg)>0.
    elif isinstance(arg,int):
       return not Lsup(arg)>0.
    elif isinstance(arg,Symbol):
       return False
    else:
       return False

def matchType(arg0=0.,arg1=0.):
    """
    converting arg0 and arg1 both to the same type L{numarray.NumArray} or L{escript.Data} or, if one of arg0 or arg1 is of type L{Symbol}, the other one to be of type L{numarray.NumArray} or L{escript.Data}.

    @param arg0: first argument 
    @type arg0: L{numarray.NumArray},L{escript.Data},C{float}, C{int}, C{Symbol}
    @param arg1: second argument 
    @type arg1: L{numarray.NumArray},L{escript.Data},C{float}, C{int}, C{Symbol}
    @return: a tuple representing arg0 and arg1 with the same type or with one of them being a L{Symbol} 
    @rtype: C{tuple} of two L{numarray.NumArray}, two L{escript.Data}, a C{Symbol} and one of the types L{numarray.NumArray} or L{escript.Data}.  
    @raise TypeError: if type of arg0 or arg1 cannot be processed
    """
    if isinstance(arg0,numarray.NumArray):
       if isinstance(arg1,numarray.NumArray):
          pass
       elif isinstance(arg1,escript.Data):
          arg0=escript.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg1=numarray.array(arg1)
       elif isinstance(arg1,int):
          arg1=numarray.array(float(arg1))
       elif isinstance(arg1,Symbol):
          pass
       else:
          raise TypeError,"function: Unknown type of second argument."    
    elif isinstance(arg0,escript.Data):
       if isinstance(arg1,numarray.NumArray):
          arg1=escript.Data(arg1,arg0.getFunctionSpace())
       elif isinstance(arg1,escript.Data):
          pass
       elif isinstance(arg1,float):
          arg1=escript.Data(arg1,(),arg0.getFunctionSpace())
       elif isinstance(arg1,int):
          arg1=escript.Data(float(arg1),(),arg0.getFunctionSpace())
       elif isinstance(arg1,Symbol):
          pass
       else:
          raise TypeError,"function: Unknown type of second argument."    
    elif isinstance(arg0,Symbol):
       if isinstance(arg1,numarray.NumArray):
          pass
       elif isinstance(arg1,escript.Data):
          pass
       elif isinstance(arg1,float):
          arg1=numarray.array(arg1)
       elif isinstance(arg1,int):
          arg1=numarray.array(float(arg1))
       elif isinstance(arg1,Symbol):
          pass
       else:
          raise TypeError,"function: Unknown type of second argument."    
    elif isinstance(arg0,float):
       if isinstance(arg1,numarray.NumArray):
          arg0=numarray.array(arg0)
       elif isinstance(arg1,escript.Data):
          arg0=escript.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numarray.array(arg0)
          arg1=numarray.array(arg1)
       elif isinstance(arg1,int):
          arg0=numarray.array(arg0)
          arg1=numarray.array(float(arg1))
       elif isinstance(arg1,Symbol):
          arg0=numarray.array(arg0)
       else:
          raise TypeError,"function: Unknown type of second argument."    
    elif isinstance(arg0,int):
       if isinstance(arg1,numarray.NumArray):
          arg0=numarray.array(float(arg0))
       elif isinstance(arg1,escript.Data):
          arg0=escript.Data(float(arg0),arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numarray.array(float(arg0))
          arg1=numarray.array(arg1)
       elif isinstance(arg1,int):
          arg0=numarray.array(float(arg0))
          arg1=numarray.array(float(arg1))
       elif isinstance(arg1,Symbol):
          arg0=numarray.array(float(arg0))
       else:
          raise TypeError,"function: Unknown type of second argument."    
    else:
      raise TypeError,"function: Unknown type of first argument."    

    return arg0,arg1

def matchShape(arg0,arg1):
    """
    

    If shape is not given the shape "largest" shape of args is used.

    @param args: a given ob
    @type arg: typically L{numarray.NumArray},L{escript.Data},C{float}, C{int}
    @return: True if the argument is identical to zero.
    @rtype: C{list} of C{int}
    """
    sh=commonShape(arg0,arg1)
    sh0=pokeShape(arg0)
    sh1=pokeShape(arg1)
    if len(sh0)<len(sh):
       return outer(arg0,numarray.ones(sh[len(sh0):],numarray.Float)),arg1
    elif len(sh1)<len(sh):
       return arg0,outer(arg1,numarray.ones(sh[len(sh1):],numarray.Float))
    else: 
       return arg0,arg1
#=========================================================
#   symbolic tool box starts here:
#=========================================================
class Symbol(object):
   """
   Symbol class. 

   Symbol class objects provide the same functionality as L{numarray.NumArray} and L{escript.Data} objects 
   but they do not have a value and therefore cannot be plotted or visualize. The main purpose is the possibilty 
   calculate derivatives with respect to other Symbols used to define a Symbol.
 
   """
   def __init__(self,shape=(),args=[],dim=None):
       """
       Creates an instance of a symbol of a given shape. The symbol may depending on a list of arguments args which may be
       symbols or any other object.

       @param arg: the arguments of the symbol.
       @type arg: C{list}
       @param shape: the shape 
       @type shape: C{tuple} of C{int}
       @param dim: spatial dimension of the symbol. If dim=C{None} the spatial dimension is undefined.  
       @type dim: C{None} or C{int}  

       """
       if len(shape)>4:
           raise ValueError,"Symbol supports only tensors up to order 4"
       self.__args=args
       self.__shape=shape
       self.__dim=dim

   def getArgument(self,i=None):
       """
       returns the i-th argument of the symbol

       @param i: index of the argument requested. 
       @type i: C{int} or C{None}
       @raise IndexError: if the requested index does not exist
       @return: the vlaue of the i-th argument or i is not specified the list of all arguments.
       @rtype: a single object or a list of objects 
       """
       if i==None:
          return self.__args
       else:
          if i<0 or i>=len(self.__args):
             raise IndexError,"there are only %s arguments"%len(self.__args)
          return self.__args[i]

   def getRank(self):
       """
       the rank of the symbol 

       @return: the rank of the symbol. This is length of the shape
       @rtype: C{int}
       """
       return len(self.getShape())

   def getShape(self):
       """
       the shape of the symbol.

       @return : the shape of the symbol. 
       @rtype: C{tuple} of C{int}
       """
       return self.__shape

   def getDim(self):
       """
       the spatial dimension

       @return : the spatial dimension
       @rtype: C{int} if the dimension is defined. Otherwise C{None} is returned.
       """
       return self.__dim
   
   def __str__(self):
       """
       a string representation of the symbol.
       @return: a string representation of the object
       @rtype: C{str}
       """
       args=[]
       for arg in self.getArgument():
          args.append(str(arg))
       try:
           out=self.getMyCode(args,format="str")
       except NotImplementedError:
           out="<Symbol %s>"%id(self)
       return out
       
   def getSubstitutedArguments(self,argvals):
       """
       substitutes symbols in the arguments of this object and returns the result as a list.

       @param argvals: L{Symbols} and their substitutes. The L{Symbol} u in the expression defining this object is replaced by argvals[u]. 
       @type argvals: C{dict} with keywords of type L{Symbol}.
       @rtype: C{list} of objects
       @return: list of the object assigned to the arguments through substitution or for the arguments which are not L{Symbols} the value assigned to the argument at instantiation.
       """
       out=[]
       for a in self.getArgument():
          if isinstance(a,Symbol):
             out.append(a.substitute(argvals))
          else:
             out.append(a)
       return out

   def getDifferentiatedArguments(self,arg):
       """
       applifies differentials to the arguments of this object and returns the result as a list.

       @param arg: the derivative is calculated with respect to arg
       @type arg: typically L{escript.Symbol} but can also be C{float}, L{escript.Data}, L{numarray.NumArray} depending the involved functions and data.
       @rtype: C{list} of objects
       @return: list of object obtained by calculating the derivatives of the argumenst with respct to arg
       """
       out=[]
       for a in self.getArgument():
          if isinstance(a,Symbol):
             out.append(a.substitute(argvals))
          else:
              s=pokeShape(s)+arg.getShape()
              if len(s)>0:
                 out.append(numarray.zeros(s),numarray.Float)
              else:
                 out.append(a)
       return out

   def isAppropriateValue(self,arg):
      """
      checks if the given argument arg can be used as a substitution of this object. The method checks 
      the shape of arg and, if the spatial dimension is defined, the spatial dimension of arg.    

      @param arg: a given object 
      @type arg: L{numarray.NumArray},L{escript.Data},C{float}, C{int}, C{Symbol}
      @return: True if arg is a suitbale object to be used for substitution. Otherwise False is returned.
      @rtype: C{bool}
      """
      if isinstance(arg,numarray.NumArray):
          return arg.shape==self.getShape()
      elif isinstance(arg,escript.Data):
          if self.getDim()==None:
              return arg.getShape()==self.getShape()
          elif self.getDim()==arg.getFunctionSpace().getDim():
              return arg.getShape()==self.getShape()
          else:
              return False
      elif isinstance(arg,Symbol):
          if self.getDim()==None:
              return arg.getShape()==self.getShape()
          elif self.getDim()==arg.getDim():
              return arg.getShape()==self.getShape()
          else:
              return False
      elif isinstance(arg,float):
          return ()==self.getShape()
      elif isinstance(arg,int):
          return ()==self.getShape()
      else:
         return False

   def getMyCode(self,argstrs,format="escript"):
       """
       returns a program code that can be used to evaluate the symbol.

       @param argstrs: gives for each argument a string representing the argument for the evaluation. 
       @type argstrs: C{list} of C{str}.
       @param format: specifies the format to be used. At the moment only "escript", "str" and "text" are supported.
       @type format: C{str}
       @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
       @rtype: C{str}
       @raise NotImplementedError: if no implementation for the given format is available
       @note: This method has to be overwritten by subclasses.
       """
       raise NotImplementedError,"no code for %s representation available"%format

   def substitute(self,argvals):
      """  
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object. 

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @note: this method has to be overwritten by a particular L{Symbol}
      @raise NotImplementedError: if no implementation for the given format is available
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"Symbol: new value is not appropriate."
      else:
         raise NotImplementedError,"no substitution in %s avialable"%str(self)

   def diff(self,arg):
       """
       returns the derivative of the symbol with respect to L{Symbol} arg

       @param arg: the derivative is calculated with respect to arg
       @type arg: typically L{escript.Symbol} but can also be C{float}, L{escript.Data}, L{numarray.NumArray} depending the involved functions and data.
       @return: derivative with respect to C{arg}
       @rtype: typically L{escript.Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
       @note: this method is overwritten by a particular L{Symbol}
       """
       if arg==self:
          return identity(self.getShape())
       else:
          s=self.getShape()+arg.getShape()
          if len(s)>0:
             return numarray.zeros(s,numarray.Float)
          else:
             return 0.

   def __neg__(self):
       """
       returns -self.

       @return:  a S{Symbol} representing the negative of the object
       @rtype: L{DependendSymbol}
       """
       return self*(-1.)

   def __pos__(self):
       """
       returns +self.

       @return:  a S{Symbol} representing the positive of the object
       @rtype: L{DependendSymbol}
       """
       return self*(1.)

   def __abs__(self):
       """
       returns a S{Symbol} representing the absolute value of the object.
       """
       return Abs_Symbol(self)

   def __add__(self,other):
       """
       add another object to this object

       @param other: object to be added to this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return:  a S{Symbol} representing the sum of this object and C{other}
       @rtype: L{DependendSymbol}
       """
       return add(self,other)

   def __radd__(self,other):
       """
       add this object to another object

       @param other: object this object is added to
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the sum of C{other} and this object object
       @rtype: L{DependendSymbol}
       """
       return add(other,self)

   def __sub__(self,other):
       """
       subtracts another object from this object

       @param other: object to be subtracted from this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the difference of C{other} and this object 
       @rtype: L{DependendSymbol}
       """
       return add(self,-other)

   def __rsub__(self,other):
       """
       subtracts this object from another object

       @param other: object this object is been subtracted from
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the difference of this object and C{other}.
       @rtype: L{DependendSymbol}
       """
       return add(-self,other)

   def __mul__(self,other):
       """
       multiplies this object with other object

       @param other: object to be mutiplied by this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the product of the object and C{other}.
       @rtype: L{DependendSymbol} or 0 if other is identical to zero.
       """
       return mult(self,other)

   def __rmul__(self,other):
       """
       multiplies this object with other object

       @param other: object this object is multiplied with
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the product of C{other} and the object.
       @rtype: L{DependendSymbol} or 0 if other is identical to zero.
       """
       return mult(other,self)

   def __div__(self,other):
       """
       divides this object by other object

       @param other: object dividing this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the quotient of this object and C{other}
       @rtype: L{DependendSymbol} 
       """
       return quotient(self,other)

   def __rdiv__(self,other):
       """
       divides this object by other object

       @param other: object dividing this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the quotient of C{other} and this object
       @rtype: L{DependendSymbol} or 0 if C{other} is identical to zero.
       """
       return quotient(other,self)

   def __pow__(self,other):
       """
       raises this object to the power of other

       @param other: exponent
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the power of this object to C{other}
       @rtype: L{DependendSymbol} or 1 if C{other} is identical to zero.
       """
       return power(self,other)

   def __rpow__(self,other):
       """
       raises an object to the power of this object

       @param other: basis 
       @type other: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @return: a S{Symbol} representing the power of C{other} to this object
       @rtype: L{DependendSymbol} or 0 if C{other} is identical to zero.
       """
       return power(other,self)

class DependendSymbol(Symbol):
   """
   DependendSymbol extents L{Symbol} by modifying the == operator to allow two instances to be equal. 
   Two DependendSymbol are equal if they have the same shape, the same arguments and one of them has an unspecified spatial dimension or the spatial dimension is identical  
   
   Example: 
   
   u1=Symbol(shape=(3,4),dim=2,args=[4.])
   u2=Symbol(shape=(3,4),dim=2,args=[4.])
   print u1==u2
   False
   
      but   

   u1=DependendSymbol(shape=(3,4),dim=2,args=[4.])
   u2=DependendSymbol(shape=(3,4),dim=2,args=[4.])
   u3=DependendSymbol(shape=(2,),dim=2,args=[4.])   
   print u1==u2, u1==u3
   True False

   @note: DependendSymbol should be used as return value of functions with L{Symbol} arguments. This will allow the optimizer to remove redundant function calls. 
   """
   def __eq__(self,other):
      """
      checks if other equals self

      @param other: any object
      @return: True if other has the same class like self, and the shape, the spatial diemsnion and the arguments are equal.
      @rtype: C{bool}
      """
      if isinstance(other,DependendSymbol):
         if self.__class__==other.__class__:  
            if self.getShape()==other.getShape():
              if self.getArgument()==other.getArgument():
                 if self.getDim()==None or other.getDim()==None or self.getDim()==other.getDim():
                    return True
      return False

   def __ne__(self,other):
      """
      checks if other equals self

      @param other: any object
      @return: Flase if other has the same class like self, and the shape, the spatial diemsnion and the arguments are equal.
      @rtype: C{bool}
      """
      return not self==other
#=========================================================
#  Unary operations prserving the shape
#========================================================
def log10(arg):
   """
   returns base-10 logarithm of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.log10(arg)
   elif isinstance(arg,escript.Data):
      return arg._log10()
   elif isinstance(arg,float):
      return math.log10(arg)
   elif isinstance(arg,int):
      return math.log10(float(arg))
   elif isinstance(arg,Symbol):
      return log(arg)/log(10.)
   else:
      raise TypeError,"log10: Unknown argument type."

def wherePositive(arg):
   """
   returns mask of positive values of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.greater(arg,numarray.zeros(arg.shape,numarray.Float))*1.
      if isinstance(out,float): out=numarray.array(out)
      return out
   elif isinstance(arg,escript.Data):
      return arg._wherePositive()
   elif isinstance(arg,float):
      if arg>0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if arg>0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return WherePositive_Symbol(arg)
   else:
      raise TypeError,"wherePositive: Unknown argument type."

class WherePositive_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the mask of positive values function
   """
   def __init__(self,arg):
      """
      initialization of wherePositive L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "wherePositive(%s)"%argstrs
      else:
         raise NotImplementedError,"WherePositive_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return wherePositive(arg)

def whereNegative(arg):
   """
   returns mask of positive values of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.less(arg,numarray.zeros(arg.shape,numarray.Float))*1.
      if isinstance(out,float): out=numarray.array(out)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNegative()
   elif isinstance(arg,float):
      if arg<0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if arg<0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return WhereNegative_Symbol(arg)
   else:
      raise TypeError,"whereNegative: Unknown argument type."

class WhereNegative_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the mask of positive values function
   """
   def __init__(self,arg):
      """
      initialization of whereNegative L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "whereNegative(%s)"%argstrs
      else:
         raise NotImplementedError,"WhereNegative_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return whereNegative(arg)

def whereNonNegative(arg):
   """
   returns mask of non-negative values of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.greater_equal(arg,numarray.zeros(arg.shape,numarray.Float))*1.
      if isinstance(out,float): out=numarray.array(out)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNonNegative()
   elif isinstance(arg,float):
      if arg<0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,int):
      if arg<0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,Symbol):
      return 1.-whereNegative(arg)
   else:
      raise TypeError,"whereNonNegative: Unknown argument type."

def whereNonPositive(arg):
   """
   returns mask of non-positive values of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.less_equal(arg,numarray.zeros(arg.shape,numarray.Float))*1.
      if isinstance(out,float): out=numarray.array(out)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNonPositive()
   elif isinstance(arg,float):
      if arg>0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,int):
      if arg>0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,Symbol):
      return 1.-wherePositive(arg)
   else:
      raise TypeError,"whereNonPositive: Unknown argument type."

def whereZero(arg,tol=0.):
   """
   returns mask of zero entries of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @param tol: tolerance. values with absolute value less then tol are accepted as zero.
   @type tol: C{float}
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.less_equal(abs(arg)-tol,numarray.zeros(arg.shape,numarray.Float))*1.
      if isinstance(out,float): out=numarray.array(out)
      return out
   elif isinstance(arg,escript.Data):
      if tol>0.:
         return whereNegative(abs(arg)-tol)
      else:
         return arg._whereZero()
   elif isinstance(arg,float):
      if abs(arg)<=tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if abs(float(arg))<=tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return WhereZero_Symbol(arg,tol)
   else:
      raise TypeError,"whereZero: Unknown argument type."

class WhereZero_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the mask of zero entries function
   """
   def __init__(self,arg,tol=0.):
      """
      initialization of whereZero L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg,tol],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "whereZero(%s,tol=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"WhereZero_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return whereZero(arg[0],arg[1])

def whereNonZero(arg,tol=0.):
   """
   returns mask of values different from zero of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.greater(abs(arg)-tol,numarray.zeros(arg.shape,numarray.Float))*1.
      if isinstance(out,float): out=numarray.array(out)
      return out
   elif isinstance(arg,escript.Data):
      if tol>0.:
         return 1.-whereZero(arg,tol)
      else:
         return arg._whereNonZero()
   elif isinstance(arg,float):
      if abs(arg)>tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if abs(float(arg))>tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return 1.-whereZero(arg,tol)
   else:
      raise TypeError,"whereNonZero: Unknown argument type."

def sin(arg):
   """
   returns sine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.sin(arg)
   elif isinstance(arg,escript.Data):
      return arg._sin()
   elif isinstance(arg,float):
      return math.sin(arg)
   elif isinstance(arg,int):
      return math.sin(arg)
   elif isinstance(arg,Symbol):
      return Sin_Symbol(arg)
   else:
      raise TypeError,"sin: Unknown argument type."

class Sin_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the sine function
   """
   def __init__(self,arg):
      """
      initialization of sin L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "sin(%s)"%argstrs
      else:
         raise NotImplementedError,"Sin_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return sin(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(cos(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def cos(arg):
   """
   returns cosine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.cos(arg)
   elif isinstance(arg,escript.Data):
      return arg._cos()
   elif isinstance(arg,float):
      return math.cos(arg)
   elif isinstance(arg,int):
      return math.cos(arg)
   elif isinstance(arg,Symbol):
      return Cos_Symbol(arg)
   else:
      raise TypeError,"cos: Unknown argument type."

class Cos_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the cosine function
   """
   def __init__(self,arg):
      """
      initialization of cos L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "cos(%s)"%argstrs
      else:
         raise NotImplementedError,"Cos_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return cos(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(-sin(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def tan(arg):
   """
   returns tangent of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.tan(arg)
   elif isinstance(arg,escript.Data):
      return arg._tan()
   elif isinstance(arg,float):
      return math.tan(arg)
   elif isinstance(arg,int):
      return math.tan(arg)
   elif isinstance(arg,Symbol):
      return Tan_Symbol(arg)
   else:
      raise TypeError,"tan: Unknown argument type."

class Tan_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the tangent function
   """
   def __init__(self,arg):
      """
      initialization of tan L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "tan(%s)"%argstrs
      else:
         raise NotImplementedError,"Tan_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return tan(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./cos(myarg)**2,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def asin(arg):
   """
   returns inverse sine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arcsin(arg)
   elif isinstance(arg,escript.Data):
      return arg._asin()
   elif isinstance(arg,float):
      return math.asin(arg)
   elif isinstance(arg,int):
      return math.asin(arg)
   elif isinstance(arg,Symbol):
      return Asin_Symbol(arg)
   else:
      raise TypeError,"asin: Unknown argument type."

class Asin_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse sine function
   """
   def __init__(self,arg):
      """
      initialization of asin L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "asin(%s)"%argstrs
      else:
         raise NotImplementedError,"Asin_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return asin(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./sqrt(1.-myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def acos(arg):
   """
   returns inverse cosine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arccos(arg)
   elif isinstance(arg,escript.Data):
      return arg._acos()
   elif isinstance(arg,float):
      return math.acos(arg)
   elif isinstance(arg,int):
      return math.acos(arg)
   elif isinstance(arg,Symbol):
      return Acos_Symbol(arg)
   else:
      raise TypeError,"acos: Unknown argument type."

class Acos_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse cosine function
   """
   def __init__(self,arg):
      """
      initialization of acos L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "acos(%s)"%argstrs
      else:
         raise NotImplementedError,"Acos_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return acos(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(-1./sqrt(1.-myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def atan(arg):
   """
   returns inverse tangent of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arctan(arg)
   elif isinstance(arg,escript.Data):
      return arg._atan()
   elif isinstance(arg,float):
      return math.atan(arg)
   elif isinstance(arg,int):
      return math.atan(arg)
   elif isinstance(arg,Symbol):
      return Atan_Symbol(arg)
   else:
      raise TypeError,"atan: Unknown argument type."

class Atan_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse tangent function
   """
   def __init__(self,arg):
      """
      initialization of atan L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "atan(%s)"%argstrs
      else:
         raise NotImplementedError,"Atan_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return atan(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./(1+myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def sinh(arg):
   """
   returns hyperbolic sine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.sinh(arg)
   elif isinstance(arg,escript.Data):
      return arg._sinh()
   elif isinstance(arg,float):
      return math.sinh(arg)
   elif isinstance(arg,int):
      return math.sinh(arg)
   elif isinstance(arg,Symbol):
      return Sinh_Symbol(arg)
   else:
      raise TypeError,"sinh: Unknown argument type."

class Sinh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the hyperbolic sine function
   """
   def __init__(self,arg):
      """
      initialization of sinh L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "sinh(%s)"%argstrs
      else:
         raise NotImplementedError,"Sinh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return sinh(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(cosh(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def cosh(arg):
   """
   returns hyperbolic cosine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.cosh(arg)
   elif isinstance(arg,escript.Data):
      return arg._cosh()
   elif isinstance(arg,float):
      return math.cosh(arg)
   elif isinstance(arg,int):
      return math.cosh(arg)
   elif isinstance(arg,Symbol):
      return Cosh_Symbol(arg)
   else:
      raise TypeError,"cosh: Unknown argument type."

class Cosh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the hyperbolic cosine function
   """
   def __init__(self,arg):
      """
      initialization of cosh L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "cosh(%s)"%argstrs
      else:
         raise NotImplementedError,"Cosh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return cosh(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(sinh(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def tanh(arg):
   """
   returns hyperbolic tangent of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.tanh(arg)
   elif isinstance(arg,escript.Data):
      return arg._tanh()
   elif isinstance(arg,float):
      return math.tanh(arg)
   elif isinstance(arg,int):
      return math.tanh(arg)
   elif isinstance(arg,Symbol):
      return Tanh_Symbol(arg)
   else:
      raise TypeError,"tanh: Unknown argument type."

class Tanh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the hyperbolic tangent function
   """
   def __init__(self,arg):
      """
      initialization of tanh L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "tanh(%s)"%argstrs
      else:
         raise NotImplementedError,"Tanh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return tanh(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./cosh(myarg)**2,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def asinh(arg):
   """
   returns inverse hyperbolic sine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arcsinh(arg)
   elif isinstance(arg,escript.Data):
      return arg._asinh()
   elif isinstance(arg,float):
      return numarray.arcsinh(arg)
   elif isinstance(arg,int):
      return numarray.arcsinh(float(arg))
   elif isinstance(arg,Symbol):
      return Asinh_Symbol(arg)
   else:
      raise TypeError,"asinh: Unknown argument type."

class Asinh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse hyperbolic sine function
   """
   def __init__(self,arg):
      """
      initialization of asinh L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "asinh(%s)"%argstrs
      else:
         raise NotImplementedError,"Asinh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return asinh(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./sqrt(myarg**2+1),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def acosh(arg):
   """
   returns inverse hyperolic cosine of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arccosh(arg)
   elif isinstance(arg,escript.Data):
      return arg._acosh()
   elif isinstance(arg,float):
      return numarray.arccosh(arg)
   elif isinstance(arg,int):
      return numarray.arccosh(float(arg))
   elif isinstance(arg,Symbol):
      return Acosh_Symbol(arg)
   else:
      raise TypeError,"acosh: Unknown argument type."

class Acosh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse hyperolic cosine function
   """
   def __init__(self,arg):
      """
      initialization of acosh L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "acosh(%s)"%argstrs
      else:
         raise NotImplementedError,"Acosh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return acosh(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./sqrt(myarg**2-1),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def atanh(arg):
   """
   returns inverse hyperbolic tangent of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arctanh(arg)
   elif isinstance(arg,escript.Data):
      return arg._atanh()
   elif isinstance(arg,float):
      return numarray.arctanh(arg)
   elif isinstance(arg,int):
      return numarray.arctanh(float(arg))
   elif isinstance(arg,Symbol):
      return Atanh_Symbol(arg)
   else:
      raise TypeError,"atanh: Unknown argument type."

class Atanh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse hyperbolic tangent function
   """
   def __init__(self,arg):
      """
      initialization of atanh L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "atanh(%s)"%argstrs
      else:
         raise NotImplementedError,"Atanh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return atanh(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./(1.-myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def exp(arg):
   """
   returns exponential of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.exp(arg)
   elif isinstance(arg,escript.Data):
      return arg._exp()
   elif isinstance(arg,float):
      return math.exp(arg)
   elif isinstance(arg,int):
      return math.exp(arg)
   elif isinstance(arg,Symbol):
      return Exp_Symbol(arg)
   else:
      raise TypeError,"exp: Unknown argument type."

class Exp_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the exponential function
   """
   def __init__(self,arg):
      """
      initialization of exp L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "exp(%s)"%argstrs
      else:
         raise NotImplementedError,"Exp_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return exp(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(self,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def sqrt(arg):
   """
   returns square root of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.sqrt(arg)
   elif isinstance(arg,escript.Data):
      return arg._sqrt()
   elif isinstance(arg,float):
      return math.sqrt(arg)
   elif isinstance(arg,int):
      return math.sqrt(arg)
   elif isinstance(arg,Symbol):
      return Sqrt_Symbol(arg)
   else:
      raise TypeError,"sqrt: Unknown argument type."

class Sqrt_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the square root function
   """
   def __init__(self,arg):
      """
      initialization of sqrt L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "sqrt(%s)"%argstrs
      else:
         raise NotImplementedError,"Sqrt_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return sqrt(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(0.5/self,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def log(arg):
   """
   returns natural logarithm of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.log(arg)
   elif isinstance(arg,escript.Data):
      return arg._log()
   elif isinstance(arg,float):
      return math.log(arg)
   elif isinstance(arg,int):
      return math.log(arg)
   elif isinstance(arg,Symbol):
      return Log_Symbol(arg)
   else:
      raise TypeError,"log: Unknown argument type."

class Log_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the natural logarithm function
   """
   def __init__(self,arg):
      """
      initialization of log L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "log(%s)"%argstrs
      else:
         raise NotImplementedError,"Log_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return log(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./arg,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def sign(arg):
   """
   returns sign of argument arg

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      return wherePositive(arg)-whereNegative(arg)
   elif isinstance(arg,escript.Data):
      return arg._sign()
   elif isinstance(arg,float):
      if arg>0:
        return 1.
      elif arg<0:
        return -1.
      else:
        return 0.
   elif isinstance(arg,int):
      if float(arg)>0:
        return 1.
      elif float(arg)<0:
        return -1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return wherePositive(arg)-whereNegative(arg)
   else:
      raise TypeError,"sign: Unknown argument type."

class Abs_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the absolute value function
   """
   def __init__(self,arg):
      """
      initialization of abs L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "abs(%s)"%argstrs
      else:
         raise NotImplementedError,"Abs_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return abs(arg)

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(sign(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def minval(arg):
   """
   returns minimum value over all components of arg at each data point

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      if arg.rank==0:
         return float(arg)
      else:
         return arg.min()
   elif isinstance(arg,escript.Data):
      return arg._minval()
   elif isinstance(arg,float):
      return arg
   elif isinstance(arg,int):
      return float(arg)
   elif isinstance(arg,Symbol):
      return Minval_Symbol(arg)
   else:
      raise TypeError,"minval: Unknown argument type."

class Minval_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the minimum value function
   """
   def __init__(self,arg):
      """
      initialization of minimum value L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "minval(%s)"%argstrs
      else:
         raise NotImplementedError,"Minval_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return minval(arg)

def maxval(arg):
   """
   returns maximum value over all components of arg at each data point

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol} depending on the type of arg.
   @raises TypeError: if the type of the argument is not expected.
   """
   if isinstance(arg,numarray.NumArray):
      if arg.rank==0:
         return float(arg)
      else:
         return arg.max()
   elif isinstance(arg,escript.Data):
      return arg._maxval()
   elif isinstance(arg,float):
      return arg
   elif isinstance(arg,int):
      return float(arg)
   elif isinstance(arg,Symbol):
      return Maxval_Symbol(arg)
   else:
      raise TypeError,"maxval: Unknown argument type."

class Maxval_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the maximum value function
   """
   def __init__(self,arg):
      """
      initialization of maximum value L{Symbol} with argument arg
      @param arg: argument of function
      @type arg: typically L{Symbol}.
      """
      DependendSymbol.__init__(self,args=[arg],shape=(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{str} or a C{list} of length 1 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript" ,"text" and "str" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "maxval(%s)"%argstrs
      else:
         raise NotImplementedError,"Maxval_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return maxval(arg)

def length(arg):
   """
   returns length/Euclidean norm of argument arg at each data point

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, L{numarray.NumArray}.
   @rtype:C{float}, L{escript.Data}, L{Symbol} depending on the type of arg.
   """
   return sqrt(inner(arg,arg))

#=======================================================
#  Binary operations:
#=======================================================
def add(arg0,arg1):
       """
       adds arg0 and arg1 together. 

       @param arg0: first term 
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: second term
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @return: the some of arg0 and arg1
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @note: The shape of both arguments is matched according to the rules in used in L{matchShape}
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]):
          return args[1]
       elif testForZero(args[1]):
          return args[0]
       else:
          if isinstance(args[0],Symbol) or isinstance(args[1],Symbol) :
              return Add_Symbol(args[0],args[1])
          elif isinstance(args[0],numarray.NumArray):
              return args[1]+args[0]
          else:
              return args[0]+args[1]

class Add_Symbol(DependendSymbol):
   """
   Symbol representing the sum of two arguments.
   """
   def __init__(self,arg0,arg1):
       """
       initialization of the L{Symbol} representing the sum of two arguments 

       @param arg0: first term in the sum
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: second term in the sum
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @raise ValueError: if both arguments do not have the same shape.
       @note: if both arguments have a spatial dimension, they must equal.
       """
       sh0=pokeShape(arg0)
       sh1=pokeShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Add_Symbol: shape of arguments must match"
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{list} of length 2 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if format=="str" or format=="text":
         return "(%s)+(%s)"%(argstrs[0],argstrs[1])
      elif format=="escript":
         return "add(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return add(args[0],args[1])

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         dargs=self.getDifferentiatedArguments(arg)
         return add(dargs[0],dargs[1])

def mult(arg0,arg1):
       """
       product of arg0 and arg1

       @param arg0: first term 
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: second term
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @return: the some of arg0 and arg1
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @note: The shape of both arguments is matched according to the rules in used in L{matchShape}
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]) or testForZero(args[1]):
          return numarray.zeros(pokeShape(args[0]),numarray.Float)
       else:
          if isinstance(args[0],Symbol) or isinstance(args[1],Symbol) :
              return Mult_Symbol(args[0],args[1])
          elif isinstance(args[0],numarray.NumArray):
              return args[1]*args[0]
          else:
              return args[0]*args[1]

class Mult_Symbol(DependendSymbol):
   """
   Symbol representing the product of two arguments.
   """
   def __init__(self,arg0,arg1):
       """
       initialization of the L{Symbol} representing the product of two arguments 

       @param arg0: first factor
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: second factor
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @raise ValueError: if both arguments do not have the same shape.
       @note: if both arguments have a spatial dimension, they must equal.
       """
       sh0=pokeShape(arg0)
       sh1=pokeShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Mult_Symbol: shape of arguments must match"
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{list} of length 2 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if format=="str" or format=="text":
         return "(%s)*(%s)"%(argstrs[0],argstrs[1])
      elif format=="escript":
         return "mult(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return mult(args[0],args[1])

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myargs=self.getArgument()
         dargs=self.getDifferentiatedArguments(arg)
         return add(mult(myargs[0],dargs[1]),mult(myargs[1],dargs[0]))

def quotient(arg0,arg1):
       """
       quotient of arg0 and arg1

       @param arg0: numerator
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: denominator
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @return: the some of arg0 and arg1
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @note: The shape of both arguments is matched according to the rules in used in L{matchShape}
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]):
          return numarray.zeros(pokeShape(args[0]),numarray.Float)
       elif isinstance(args[0],Symbol):
          if isinstance(args[1],Symbol):
             return Quotient_Symbol(args[0],args[1])
          else:
             return mult(args[0],1./args[1])
       else:
          if isinstance(args[1],Symbol):
             return Quotient_Symbol(args[0],args[1])
          elif isinstance(args[0],numarray.NumArray) and not isinstance(args[1],numarray.NumArray):
             return 1./args[1]*args[0]
          else:
             return args[0]/args[1]

class Quotient_Symbol(DependendSymbol):
   """
   Symbol representing the quotient of two arguments.
   """
   def __init__(self,arg0,arg1):
       """
       initialization of L{Symbol} representing the quotient of two arguments 

       @param arg0: numerator
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: denominator
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @raise ValueError: if both arguments do not have the same shape.
       @note: if both arguments have a spatial dimension, they must equal.
       """
       sh0=pokeShape(arg0)
       sh1=pokeShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Quotient_Symbol: shape of arguments must match"
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{list} of length 2 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if format=="str" or format=="text":
         return "(%s)/(%s)"%(argstrs[0],argstrs[1])
      if format=="escript":
         return "quotient(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return quotient(args[0],args[1])

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myargs=self.getArgument()
         dargs=self.getDifferentiatedArguments(arg)
         return quotient(add(mult(myargs[1],dargs[0]),mult(-myargs[0],dargs[1])),myargs[1]*myargs[1])


def power(arg0,arg1):
       """
       raises arg0 to the power of arg1

       @param arg0: basis
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: exponent
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @return: of arg0 and arg1
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data}, L{numarray.NumArray}.
       @note: The shape of both arguments is matched according to the rules in used in L{matchShape}
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]):
          return numarray.zeros(args[0],numarray.Float)
       elif testForZero(args[1]):
          return numarray.ones(args[0],numarray.Float)
       elif isinstance(args[0],Symbol) or isinstance(args[1],Symbol):
          return Power_Symbol(args[0],args[1])
       elif isinstance(args[0],numarray.NumArray) and not isinstance(args[1],numarray.NumArray):
          return exp(args[1]*log(args[0]))
       else:
           return args[0]**args[1]

class Power_Symbol(DependendSymbol):
   """
   Symbol representing the first argument to the power of the second argument.
   """
   def __init__(self,arg0,arg1):
       """
       initialization of the L{Symbol} representing rasing the first argument to the power of the second.

       @param arg0: basis
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: exponent
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @raise ValueError: if both arguments do not have the same shape.
       @note: if both arguments have a spatial dimension, they must equal.
       """
       sh0=pokeShape(arg0)
       sh1=pokeShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Power_Symbol: shape of arguments must match"
       d0=pokeDim(arg0)
       d1=pokeDim(arg1)
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{list} of length 2 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if format=="escript" or format=="str" or format=="text":
         return "(%s)**(%s)"%(argstrs[0],argstrs[1])
      elif format=="escript":
         return "power(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return power(args[0],args[1])

   def diff(self,arg):
      """
      differential of this object

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float}, L{escript.Data}, L{numarray.NumArray}  are possible.
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myargs=self.getArgument()
         dargs=self.getDifferentiatedArguments(arg)
         return mult(self,add(mult(log(myargs[0]),dargs[1]),mult(quotient(myargs[1],myargs[0]),dargs[0])))

def maximum(*args):
    """
    the maximum over arguments args
 
    @param args: arguments
    @type args: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or C{float}
    @return: is on object which gives at each entry the maximum of the coresponding values all args
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or C{float} depending on the input
    """
    out=None
    for a in args:
       if out==None:
          out=a
       else:
          diff=add(a,-out)
          out=add(out,mult(wherePositive(diff),diff))
    return out
   
def minimum(*args):
    """
    the minimum over arguments args
 
    @param args: arguments
    @type args: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or C{float}
    @return: is on object which gives at each entry the minimum of the coresponding values all args
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or C{float} depending on the input
    """
    out=None
    for a in args:
       if out==None:
          out=a
       else:
          diff=add(a,-out)
          out=add(out,mult(whereNegative(diff),diff))
    return out

def clip(arg,minval=0.,maxval=1.):
    """
    cuts the values of arg between minval and maxval
 
    @param arg: argument
    @type arg: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or C{float}
    @param minval: lower range 
    @type arg: C{float}
    @param maxval: upper range 
    @type arg: C{float}
    @return: is on object with all its value between minval and maxval. value of the argument that greater then minval and
             less then maxval are unchanged.
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or C{float} depending on the input
    @raise ValueError: if minval>maxval
    """
    if minval>maxval:
       raise ValueError,"minval = %s must be less then maxval %s"%(minval,maxval)
    return minimum(maximum(minval,arg),maxval)

   
def inner(arg0,arg1):
    """
    inner product of the two argument:
 
    out=S{Sigma}_s arg0[s]*arg1[s]

    where s runs through arg0.Shape.

    arg0 and arg1 must have the same shape.

    @param arg0: first argument 
    @type arg0: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}, C{int}
    @param arg1: second argument 
    @type arg1: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}, C{int}
    @return : the inner product of arg0 and arg1 at each data point
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float} depending on the input
    @raise ValueError: if the shapes of the arguments are not identical
    """
    sh0=pokeShape(arg0)
    sh1=pokeShape(arg1)
    if not sh0==sh1:
        raise ValueError,"inner: shape of arguments does not match"
    return generalTensorProduct(arg0,arg1,offset=len(sh0))

def matrixmult(arg0,arg1):
    """
    matrix-matrix or matrix-vector product of the two argument:
 
    out[s0]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0] 

          or 

    out[s0,s1]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0,s1] 
 
    The second dimension of arg0 and the length of arg1 must match.

    @param arg0: first argument of rank 2
    @type arg0: L{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of at least rank 1
    @type arg1: L{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the matrix-matrix or matrix-vector product of arg0 and arg1 at each data point
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the input
    @raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=pokeShape(arg0)
    sh1=pokeShape(arg1)
    if not len(sh0)==2 :
        raise ValueError,"first argument must have rank 2"
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError,"second argument must have rank 1 or 2"
    return generalTensorProduct(arg0,arg1,offset=1)

def outer(arg0,arg1):
    """
    the outer product of the two argument:
 
    out[t,s]=arg0[t]*arg1[s]

    where s runs through arg0.Shape
          t runs through arg1.Shape

    @param arg0: first argument 
    @type arg0: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}, C{int}
    @param arg1: second argument 
    @type arg1: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}, C{int}
    @return: the outer product of arg0 and arg1 at each data point
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the input
    """
    return generalTensorProduct(arg0,arg1,offset=0)


def tensormult(arg0,arg1):
    """
    the tensor product of the two argument:

    
    for arg0 of rank 2 this is
 
    out[s0]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0]  

                 or 

    out[s0,s1]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0,s1] 

   
    and for arg0 of rank 4 this is 

    out[s0,s1,s2,s3]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2,s3] 
              
                 or

    out[s0,s1,s2]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2] 
 
                 or 

    out[s0,s1]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1] 

    In the first case the the second dimension of arg0 and the length of arg1 must match and  
    in the second case the two last dimensions of arg0 must match the shape of arg1.

    @param arg0: first argument of rank 2 or 4
    @type arg0: L{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of shape greater of 1 or 2 depending on rank of arg0
    @type arg1: L{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the tensor product of arg0 and arg1 at each data point
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the input
    """
    sh0=pokeShape(arg0)
    sh1=pokeShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
       return generalTensorProduct(arg0,arg1,offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
       return generalTensorProduct(arg0,arg1,offset=2)
    else:
        raise ValueError,"tensormult: first argument must have rank 2 or 4"

def generalTensorProduct(arg0,arg1,offset=0):
    """
    generalized tensor product 

    out[s,t]=S{Sigma}_r arg0[s,r]*arg1[r,t]

    where s runs through arg0.Shape[:arg0.Rank-offset]
          r runs trough arg0.Shape[:offset]
          t runs through arg1.Shape[offset:]

    In the first case the the second dimension of arg0 and the length of arg1 must match and  
    in the second case the two last dimensions of arg0 must match the shape of arg1.

    @param arg0: first argument
    @type arg0: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}, C{int}
    @param arg1: second argument of shape greater of 1 or 2 depending on rank of arg0
    @type arg1: L{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}, C{int}
    @return: the general tensor product of arg0 and arg1 at each data point. 
    @rtype: L{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the input
    """
    if isinstance(arg0,float) and isinstance(arg1,float): return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg0 are both numarray.NumArray or escript.Data or Symbols
    if isinstance(arg0,numarray.NumArray):
       if isinstance(arg1,Symbol):
           return GeneralTensorProduct_Symbol(arg0,arg1,offset)
       else:
           if not arg0.shape[arg0.rank-offset:]==arg1.shape[:offset]:
               raise ValueError,"generalTensorProduct: dimensions of last %s components in left argument don't match the first %s components in the right argument."%(offset,offset) 
           arg0_c=arg0.copy()
           arg1_c=arg1.copy()
           sh0,sh1=arg0.shape,arg1.shape
           d0,d1,d01=1,1,1
           for i in sh0[:arg0.rank-offset]: d0*=i
           for i in sh1[offset:]: d1*=i
           for i in sh1[:offset]: d01*=i
           arg0_c.resize((d0,d01))
           arg1_c.resize((d01,d1))
           out=numarray.zeros((d0,d1),numarray.Float)
           for i0 in range(d0):
                    for i1 in range(d1):
                         out[i0,i1]=numarray.sum(arg0_c[i0,:]*arg1_c[:,i1])
           out.resize(sh0[:arg0.rank-offset]+sh1[offset:])
           return out
    elif isinstance(arg0,escript.Data):
       if isinstance(arg1,Symbol):
           return GeneralTensorProduct_Symbol(arg0,arg1,offset)
       else:
           return escript_generalTensorProduct(arg0,arg1,offset) # this calls has to be replaced by escript._generalTensorProduct(arg0,arg1,offset)
    else:       
       return GeneralTensorProduct_Symbol(arg0,arg1,offset)
                 
class GeneralTensorProduct_Symbol(DependendSymbol):
   """
   Symbol representing the quotient of two arguments.
   """
   def __init__(self,arg0,arg1,offset=0):
       """
       initialization of L{Symbol} representing the quotient of two arguments 

       @param arg0: numerator
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @param arg1: denominator
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray}.
       @raise ValueError: if both arguments do not have the same shape.
       @note: if both arguments have a spatial dimension, they must equal.
       """
       sh_arg0=pokeShape(arg0)
       sh_arg1=pokeShape(arg1)
       sh0=sh_arg0[:len(sh_arg0)-offset]
       sh01=sh_arg0[len(sh_arg0)-offset:]
       sh10=sh_arg1[:offset]
       sh1=sh_arg1[offset:]
       if not sh01==sh10:
           raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(offset,offset)
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0+sh1,args=[arg0,arg1,offset])

   def getMyCode(self,argstrs,format="escript"):
      """
      returns a program code that can be used to evaluate the symbol.

      @param argstrs: gives for each argument a string representing the argument for the evaluation.
      @type argstrs: C{list} of length 2 of C{str}.
      @param format: specifies the format to be used. At the moment only "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the expression assuming the values for the arguments are available.
      @rtype: C{str}
      @raise: NotImplementedError: if the requested format is not available
      """
      if format=="escript" or format=="str" or format=="text":
         return "generalTensorProduct(%s,%s,offset=%s)"%(argstrs[0],argstrs[1],argstrs[2])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, L{numarray.NumArray} depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted.
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return generalTensorProduct(args[0],args[1],args[2])

def escript_generalTensorProduct(arg0,arg1,offset): # this should be escript._generalTensorProduct
    "arg0 and arg1 are both Data objects but not neccesrily on the same function space. they could be identical!!!"
    # calculate the return shape:
    shape0=arg0.getShape()[:arg0.getRank()-offset]
    shape01=arg0.getShape()[arg0.getRank()-offset:]
    shape10=arg1.getShape()[:offset]
    shape1=arg1.getShape()[offset:]
    if not shape01==shape10:
        raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(offset,offset) 

    # whatr function space should be used? (this here is not good!)
    fs=(escript.Scalar(0.,arg0.getFunctionSpace())+escript.Scalar(0.,arg1.getFunctionSpace())).getFunctionspace()
    # create return value:
    out=escript.Data(0.,tuple(shape0+shape1),fs)
    # 
    s0=[[]]
    for k in shape0: 
          s=[]
          for j in s0: 
                for i in range(k): s.append(j+[slice(i,i)])
          s0=s
    s1=[[]]
    for k in shape1: 
          s=[]
          for j in s1: 
                for i in range(k): s.append(j+[slice(i,i)])
          s1=s
    s01=[[]]
    for k in shape01: 
          s=[]
          for j in s01: 
                for i in range(k): s.append(j+[slice(i,i)])
          s01=s

    for i0 in s0:
       for i1 in s1:
         s=escript.Scalar(0.,fs)
         for i01 in s01:
            s+=arg0.__getitem__(tuple(i0+i01))*arg1.__getitem__(tuple(i01+i1))
         out.__setitem__(tuple(i0+i1),s)
    return out

#=========================================================
#   some little helpers
#=========================================================
def grad(arg,where=None):
    """
    Returns the spatial gradient of arg at where.

    @param arg:   Data object representing the function which gradient 
                  to be calculated.
    @param where: FunctionSpace in which the gradient will be calculated. 
                  If not present or C{None} an appropriate default is used. 
    """
    if isinstance(arg,Symbol):
       return Grad_Symbol(arg,where)
    elif isinstance(arg,escript.Data):
       if where==None:
          return arg._grad()
       else:
          return arg._grad(where)
    else:
      raise TypeError,"grad: Unknown argument type."

def integrate(arg,where=None):
    """
    Return the integral if the function represented by Data object arg over 
    its domain.

    @param arg:   Data object representing the function which is integrated.
    @param where: FunctionSpace in which the integral is calculated. 
                  If not present or C{None} an appropriate default is used.
    """
    if isinstance(arg,numarray.NumArray):
        if checkForZero(arg):
           return arg
        else:
           raise TypeError,"integrate: cannot intergrate argument"
    elif isinstance(arg,float):
        if checkForZero(arg):
           return arg
        else:
           raise TypeError,"integrate: cannot intergrate argument"
    elif isinstance(arg,int):
        if checkForZero(arg):
           return float(arg)
        else:
           raise TypeError,"integrate: cannot intergrate argument"
    elif isinstance(arg,Symbol):
       return Integrate_Symbol(arg,where)
    elif isinstance(arg,escript.Data):
       if not where==None: arg=escript.Data(arg,where)
       if arg.getRank()==0:
          return arg._integrate()[0]
       else:
          return arg._integrate()
    else:
      raise TypeError,"integrate: Unknown argument type."

def interpolate(arg,where):
    """
    Interpolates the function into the FunctionSpace where.

    @param arg:    interpolant 
    @param where:  FunctionSpace to interpolate to
    """
    if isinstance(arg,Symbol):
       return Interpolated_Symbol(arg,where)
    else:
       return escript.Data(arg,where)

def div(arg,where=None):
    """
    Returns the divergence of arg at where.

    @param arg:   Data object representing the function which gradient to 
                  be calculated.
    @param where: FunctionSpace in which the gradient will be calculated. 
                  If not present or C{None} an appropriate default is used. 
    """
    g=grad(arg,where)
    return trace(g,axis0=g.getRank()-2,axis1=g.getRank()-1)

def jump(arg):
    """
    Returns the jump of arg across a continuity.

    @param arg:   Data object representing the function which gradient 
                  to be calculated.
    """
    d=arg.getDomain()
    return arg.interpolate(escript.FunctionOnContactOne(d))-arg.interpolate(escript.FunctionOnContactZero(d))

#=============================
#
# wrapper for various functions: if the argument has attribute the function name
# as an argument it calls the corresponding methods. Otherwise the corresponding
# numarray function is called.

# functions involving the underlying Domain:

def transpose(arg,axis=None):
    """
    Returns the transpose of the Data object arg. 

    @param arg:
    """
    if axis==None:
       r=0
       if hasattr(arg,"getRank"): r=arg.getRank()
       if hasattr(arg,"rank"): r=arg.rank
       axis=r/2
    if isinstance(arg,Symbol):
       return Transpose_Symbol(arg,axis=r)
    if isinstance(arg,escript.Data):
       # hack for transpose 
       r=arg.getRank()
       if r!=2: raise ValueError,"Tranpose only avalaible for rank 2 objects"
       s=arg.getShape()
       out=escript.Data(0.,(s[1],s[0]),arg.getFunctionSpace())
       for i in range(s[0]):
          for j in range(s[1]):
             out[j,i]=arg[i,j]
       return out
       # end hack for transpose 
       return arg.transpose(axis)
    else:
       return numarray.transpose(arg,axis=axis)

def trace(arg,axis0=0,axis1=1):
    """
    Return 

    @param arg:
    """
    if isinstance(arg,Symbol):
       s=list(arg.getShape())        
       s=tuple(s[0:axis0]+s[axis0+1:axis1]+s[axis1+1:])
       return Trace_Symbol(arg,axis0=axis0,axis1=axis1)
    elif isinstance(arg,escript.Data):
       # hack for trace 
       s=arg.getShape()
       if s[axis0]!=s[axis1]:
           raise ValueError,"illegal axis in trace"
       out=escript.Scalar(0.,arg.getFunctionSpace())
       for i in range(s[axis0]):
          out+=arg[i,i]
       return out
       # end hack for trace
    else:
       return numarray.trace(arg,axis0=axis0,axis1=axis1)


def reorderComponents(arg,index):
    """
    resorts the component of arg according to index

    """
    pass
#
# $Log: util.py,v $
# Revision 1.14.2.16  2005/10/19 06:09:57  gross
# new versions of saveVTK and saveDX which allow to write several Data objects into a single file
#
# Revision 1.14.2.15  2005/10/19 03:25:27  jgs
# numarray uses log10 for log, and log for ln - log()/ln() updated to reflect this
#
# Revision 1.14.2.14  2005/10/19 02:10:23  jgs
# fixed ln() to call Data::ln
#
# Revision 1.14.2.13  2005/09/12 03:32:14  gross
# test_visualiztion has been aded to mk
#
# Revision 1.14.2.12  2005/09/09 01:56:24  jgs
# added implementations of acos asin atan sinh cosh tanh asinh acosh atanh
# and some associated testing
#
# Revision 1.14.2.11  2005/09/08 08:28:39  gross
# some cleanup in savevtk
#
# Revision 1.14.2.10  2005/09/08 00:25:32  gross
# test for finley mesh generators added
#
# Revision 1.14.2.9  2005/09/07 10:32:05  gross
# Symbols removed from util and put into symmbols.py.
#
# Revision 1.14.2.8  2005/08/26 05:06:37  cochrane
# Corrected errors in docstrings.  Improved output formatting of docstrings.
# Other minor improvements to code and docs (eg spelling etc).
#
# Revision 1.14.2.7  2005/08/26 04:45:40  cochrane
# Fixed and tidied markup and docstrings.  Some *Symbol classes were defined
# as functions, so changed them to classes (hopefully this was the right thing
# to do).
#
# Revision 1.14.2.6  2005/08/26 04:30:13  gross
# gneric unit testing for linearPDE
#
# Revision 1.14.2.5  2005/08/24 02:02:52  gross
# jump function added
#
# Revision 1.14.2.4  2005/08/18 04:39:32  gross
# the constants have been removed from util.py as they not needed anymore. PDE related constants are accessed through LinearPDE attributes now
#
# Revision 1.14.2.3  2005/08/03 09:55:33  gross
# ContactTest is passing now./mk install!
#
# Revision 1.14.2.2  2005/08/02 03:15:14  gross
# bug inb trace fixed!
#
# Revision 1.14.2.1  2005/07/29 07:10:28  gross
# new functions in util and a new pde type in linearPDEs
#
# Revision 1.2.2.21  2005/07/28 04:19:23  gross
# new functions maximum and minimum introduced.
#
# Revision 1.2.2.20  2005/07/25 01:26:27  gross
# bug in inner fixed
#
# Revision 1.2.2.19  2005/07/21 04:01:28  jgs
# minor comment fixes
#
# Revision 1.2.2.18  2005/07/21 01:02:43  jgs
# commit ln() updates to development branch version
#
# Revision 1.12  2005/07/20 06:14:58  jgs
# added ln(data) style wrapper for data.ln() - also added corresponding
# implementation of Ln_Symbol class (not sure if this is right though)
#
# Revision 1.11  2005/07/08 04:07:35  jgs
# Merge of development branch back to main trunk on 2005-07-08
#
# Revision 1.10  2005/06/09 05:37:59  jgs
# Merge of development branch back to main trunk on 2005-06-09
#
# Revision 1.2.2.17  2005/07/07 07:28:58  gross
# some stuff added to util.py to improve functionality
#
# Revision 1.2.2.16  2005/06/30 01:53:55  gross
# a bug in coloring fixed
#
# Revision 1.2.2.15  2005/06/29 02:36:43  gross
# Symbols have been introduced and some function clarified. needs much more work
#
# Revision 1.2.2.14  2005/05/20 04:05:23  gross
# some work on a darcy flow started
#
# Revision 1.2.2.13  2005/03/16 05:17:58  matt
# Implemented unit(idx, dim) to create cartesian unit basis vectors to
# complement kronecker(dim) function.
#
# Revision 1.2.2.12  2005/03/10 08:14:37  matt
# Added non-member Linf utility function to complement Data::Linf().
#
# Revision 1.2.2.11  2005/02/17 05:53:25  gross
# some bug in saveDX fixed: in fact the bug was in
# DataC/getDataPointShape
#
# Revision 1.2.2.10  2005/01/11 04:59:36  gross
# automatic interpolation in integrate switched off
#
# Revision 1.2.2.9  2005/01/11 03:38:13  gross
# Bug in Data.integrate() fixed for the case of rank 0. The problem is not totallly resolved as the method should return a scalar rather than a numarray object in the case of rank 0. This problem is fixed by the util.integrate wrapper.
#
# Revision 1.2.2.8  2005/01/05 04:21:41  gross
# FunctionSpace checking/matchig in slicing added
#
# Revision 1.2.2.7  2004/12/29 05:29:59  gross
# AdvectivePDE successfully tested for Peclet number 1000000. there is still a problem with setValue and Data()
#
# Revision 1.2.2.6  2004/12/24 06:05:41  gross
# some changes in linearPDEs to add AdevectivePDE
#
# Revision 1.2.2.5  2004/12/17 00:06:53  gross
# mk sets ESYS_ROOT is undefined
#
# Revision 1.2.2.4  2004/12/07 03:19:51  gross
# options for GMRES and PRES20 added
#
# Revision 1.2.2.3  2004/12/06 04:55:18  gross
# function wraper extended
#
# Revision 1.2.2.2  2004/11/22 05:44:07  gross
# a few more unitary functions have been added but not implemented in Data yet
#
# Revision 1.2.2.1  2004/11/12 06:58:15  gross
# a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
#
# Revision 1.2  2004/10/27 00:23:36  jgs
# fixed minor syntax error
#
# Revision 1.1.1.1  2004/10/26 06:53:56  jgs
# initial import of project esys2
#
# Revision 1.1.2.3  2004/10/26 06:43:48  jgs
# committing Lutz's and Paul's changes to brach jgs
#
# Revision 1.1.4.1  2004/10/20 05:32:51  cochrane
# Added incomplete Doxygen comments to files, or merely put the docstrings that already exist into Doxygen form.
#
# Revision 1.1  2004/08/05 03:58:27  gross
# Bug in Assemble_NodeCoordinates fixed
#
#

# vim: expandtab shiftwidth=4:
