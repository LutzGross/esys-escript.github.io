# $Id$

## @file util.py

"""
@brief Utility functions for escript
"""

import numarray
import escript
#
#   escript constants (have to be consistent witj utilC.h 
#
UNKNOWN=-1
EPSILON=1.e-15
Pi=numarray.pi
# some solver options:
NO_REORDERING=0
MINIMUM_FILL_IN=1
NESTED_DISSECTION=2
# solver methods
DEFAULT_METHOD=0
DIRECT=1
CHOLEVSKY=2
PCG=3
CR=4
CGS=5
BICGSTAB=6
SSOR=7
ILU0=8
ILUT=9
JACOBI=10
GMRES=11
PRES20=12

METHOD_KEY="method"
SYMMETRY_KEY="symmetric"
TOLERANCE_KEY="tolerance"

# supported file formats:
VRML=1
PNG=2
JPEG=3
JPG=3
PS=4
OOGL=5
BMP=6
TIFF=7
OPENINVENTOR=8
RENDERMAN=9
PNM=10

#
# wrapper for various functions: if the argument has attribute the function name
# as an argument it calls the correspong methods. Otherwise the coresponsing numarray
# function is called.
#
# functions involving the underlying Domain:
#
def grad(arg,where=None):
    """
    @brief returns the spatial gradient of the Data object arg

    @param arg: Data object representing the function which gradient to be calculated.
    @param where: FunctionSpace in which the gradient will be. If None Function(dom) where dom is the
                  domain of the Data object arg.
    """
    if isinstance(arg,escript.Data):
       if where==None:
          return arg.grad()
       else:
          return arg.grad(where)
    else:
       return arg*0.

def integrate(arg,what=None):
    """
    @brief return the integral if the function represented by Data object arg over its domain.

    @param arg
    """
    if not what==None:
       arg2=escript.Data(arg,what)
    else:
       arg2=arg
    if arg2.getRank()==0:
        return arg2.integrate()[0]
    else:
        return arg2.integrate()

def interpolate(arg,where):
    """
    @brief interpolates the function represented by Data object arg into the FunctionSpace where.

    @param arg
    @param where
    """
    if isinstance(arg,escript.Data):
       return arg.interpolate(where)
    else:
       return arg

# functions returning Data objects:

def transpose(arg,axis=None):
    """
    @brief returns the transpose of the Data object arg. 

    @param arg
    """
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
       if axis==None: axis=arg.getRank()/2
       return arg.transpose(axis)
    else:
       if axis==None: axis=arg.rank/2
       return numarray.transpose(arg,axis=axis)

def trace(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       # hack for trace 
       r=arg.getRank()
       if r!=2: raise ValueError,"trace only avalaible for rank 2 objects"
       s=arg.getShape()
       out=escript.Scalar(0,arg.getFunctionSpace())
       for i in range(min(s)):
             out+=arg[i,i]
       return out
       # end hack for trace
       return arg.trace()
    else:
       return numarray.trace(arg)

def exp(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.exp()
    else:
       return numarray.exp(arg)

def sqrt(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.sqrt()
    else:
       return numarray.sqrt(arg)

def sin(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.sin()
    else:
       return numarray.sin(arg)

def tan(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.tan()
    else:
       return numarray.tan(arg)

def cos(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.cos()
    else:
       return numarray.cos(arg)

def maxval(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.maxval()
    elif isinstance(arg,float) or isinstance(arg,int):
       return arg
    else:
       return arg.max()

def minval(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.minval()
    elif isinstance(arg,float) or isinstance(arg,int):
       return arg
    else:
       return arg.min()

def length(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       if arg.isEmpty(): return escript.Data()
       if arg.getRank()==0:
          return abs(arg)
       elif arg.getRank()==1:
          sum=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             sum+=arg[i]**2
          return sqrt(sum)
       elif arg.getRank()==2:
          sum=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             for j in range(arg.getShape()[1]):
                sum+=arg[i,j]**2
          return sqrt(sum)
       elif arg.getRank()==3:
          sum=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             for j in range(arg.getShape()[1]):
                for k in range(arg.getShape()[2]):
                   sum+=arg[i,j,k]**2
          return sqrt(sum)
       elif arg.getRank()==4:
          sum=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             for j in range(arg.getShape()[1]):
                for k in range(arg.getShape()[2]):
                   for l in range(arg.getShape()[3]):
                      sum+=arg[i,j,k,l]**2
          return sqrt(sum)
       else:
          raise SystemError,"length is not been implemented yet"
       # return arg.length()
    else:
       return sqrt((arg**2).sum())

def deviator(arg):
    """
    @brief

    @param arg1
    """
    if isinstance(arg,escript.Data):
        shape=arg.getShape()
    else:
        shape=arg.shape
    if len(shape)!=2:
          raise ValueError,"Deviator requires rank 2 object"
    if shape[0]!=shape[1]:
          raise ValueError,"Deviator requires a square matrix"
    return arg-1./(shape[0]*1.)*trace(arg)*kronecker(shape[0])

def inner(arg1,arg2):
    """
    @brief

    @param arg1, arg2
    """
    sum=escript.Scalar(0,arg1.getFunctionSpace())
    if arg.getRank()==0:
          return arg1*arg2
    elif arg.getRank()==1:
         sum=escript.Scalar(0,arg.getFunctionSpace())
         for i in range(arg.getShape()[0]):
            sum+=arg1[i]*arg2[i]
    elif arg.getRank()==2:
        sum=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
           for j in range(arg.getShape()[1]):
              sum+=arg1[i,j]*arg2[i,j]
    elif arg.getRank()==3:
        sum=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
            for j in range(arg.getShape()[1]):
               for k in range(arg.getShape()[2]):
                  sum+=arg1[i,j,k]*arg2[i,j,k]
    elif arg.getRank()==4:
        sum=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
           for j in range(arg.getShape()[1]):
              for k in range(arg.getShape()[2]):
                 for l in range(arg.getShape()[3]):
                    sum+=arg1[i,j,k,l]*arg2[i,j,k,l]
    else:
          raise SystemError,"inner is not been implemented yet"
    return sum

def sign(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.sign()
    else:
       return numarray.greater(arg,numarray.zeros(arg.shape))-numarray.less(arg,numarray.zeros(arg.shape))

# reduction operations:

def sum(arg):
    """
    @brief

    @param arg
    """
    return arg.sum()

def sup(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.sup()
    elif isinstance(arg,float) or isinstance(arg,int):
       return arg
    else:
       return arg.max()

def inf(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.inf()
    elif isinstance(arg,float) or isinstance(arg,int):
       return arg
    else:
       return arg.min()

def L2(arg):
    """
    @brief returns the L2-norm of the 

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.L2()
    elif isinstance(arg,float) or isinstance(arg,int):
       return abs(arg)
    else:
       return numarry.sqrt(dot(arg,arg))

def Lsup(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.Lsup()
    elif isinstance(arg,float) or isinstance(arg,int):
       return abs(arg)
    else:
       return max(numarray.abs(arg))

def Linf(arg):
    """
    @brief

    @param arg
    """
    if isinstance(arg,escript.Data):
       return arg.Linf()
    elif isinstance(arg,float) or isinstance(arg,int):
       return abs(arg)
    else:
       return min(numarray.abs(arg))

def dot(arg1,arg2):
    """
    @brief

    @param arg
    """
    if isinstance(arg1,escript.Data):
       return arg1.dot(arg2)
    elif isinstance(arg1,escript.Data):
       return arg2.dot(arg1)
    else:
       return numarray.dot(arg1,arg2)

def kronecker(d):
   return numarray.identity(d)

def unit(i,d):
   """
   @brief return a unit vector of dimension d with nonzero index i
   @param d dimension
   @param i index
   """
   e = numarray.zeros((d,))
   e[i] = 1.0
   return e
