# $Id$

## @file util.py

"""
@brief Utility functions for escript

TODO for Data:

  * binary operations @:               (a@b)[:,*]=a[:]@b[:,*] if rank(a)<rank(b)
                @=+,-,*,/,**           (a@b)[:]=a[:]@b[:] if rank(a)=rank(b)
                                       (a@b)[*,:]=a[*,:]@b[:] if rank(a)>rank(b)
   
  * implementation of outer outer(a,b)[:,*]=a[:]*b[*]
  * trace: trace(arg,axis0=a0,axis1=a1)(:,&,*)=sum_i trace(:,i,&,i,*) (i are at index a0 and a1)
"""

import numarray
import escript

#
#   escript constants (have to be consistent with utilC.h )
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

#===========================================================
# a simple tool box to deal with _differentials of functions 
#===========================================================

class Symbol:
   """symbol class"""
   def __init__(self,name="symbol",shape=(),dim=3,args=[]):
       """creates an instance of a symbol of shape shape and spatial dimension dim
          The symbol may depending on a list of arguments args which
          may be symbols or other objects. name gives the name of the symbol."""
       self.__args=args
       self.__name=name
       self.__shape=shape
       if hasattr(dim,"getDim"):
           self.__dim=dim.getDim()
       else:    
           self.__dim=dim
       # 
       self.__cache_val=None
       self.__cache_argval=None

   def getArgument(self,i):
       """returns the i-th argument"""
       return self.__args[i]

   def getDim(self):
       """returns the spatial dimension of the symbol"""
       return self.__dim

   def getRank(self):
       """returns the rank of the symbol"""
       return len(self.getShape())

   def getShape(self):
       """returns the shape of the symbol"""
       return self.__shape

   def getEvaluatedArguments(self,argval):
       """returns the list of evaluated arguments by subsituting symbol u by argval[u]."""
       if argval==self.__cache_argval:
           print "%s: cached value used"%self
           return self.__cache_val
       else: 
           out=[]
           for a  in self.__args:
              if isinstance(a,Symbol):
                out.append(a.eval(argval))
              else:
                out.append(a)
           self.__cache_argval=argval
           self.__cache_val=out
           return out

   def getDifferentiatedArguments(self,arg):
       """returns the list of the arguments _differentiated by arg""" 
       out=[]
       for a in self.__args:
          if isinstance(a,Symbol):
            out.append(a.diff(arg))
          else:
            out.append(0)
       return out

   def diff(self,arg):
       """returns the _differention of self by arg."""
       if self==arg:
          out=numarray.zeros(tuple(2*list(self.getShape())),numarray.Float)
          if self.getRank()==0:
             out=1.
          elif self.getRank()==1:
              for i0 in range(self.getShape()[0]):
                 out[i0,i0]=1.   
          elif self.getRank()==2:
              for i0 in range(self.getShape()[0]):
                for i1 in range(self.getShape()[1]):
                     out[i0,i1,i0,i1]=1.   
          elif self.getRank()==3:
              for i0 in range(self.getShape()[0]):
                for i1 in range(self.getShape()[1]):
                  for i2 in range(self.getShape()[2]):
                     out[i0,i1,i2,i0,i1,i2]=1.   
          elif self.getRank()==4:
              for i0 in range(self.getShape()[0]):
                for i1 in range(self.getShape()[1]):
                  for i2 in range(self.getShape()[2]):
                    for i3 in range(self.getShape()[3]): 
                       out[i0,i1,i2,i3,i0,i1,i2,i3]=1.   
          else:
             raise ValueError,"differential support rank<5 only."
          return out
       else:
          return self._diff(arg)

   def _diff(self,arg):
       """return derivate of self with respect to arg (!=self). 
          This method is overwritten by a particular symbol"""
       return 0

   def eval(self,argval):
       """subsitutes symbol u in self by argval[u] and returns the result. If
          self is not a key of argval then self is returned."""
       if argval.has_key(self):
         return argval[self]
       else:
         return self

   def __str__(self):
       """returns a string representation of the symbol"""
       return self.__name

   def __add__(self,other):
       """adds other to symbol self. if _testForZero(other) self is returned."""
       if _testForZero(other):
          return self
       else:
          a=_matchShape([self,other])
          return Add_Symbol(a[0],a[1])

   def __radd__(self,other):
       """adds other to symbol self. if _testForZero(other) self is returned."""
       return self+other

   def __neg__(self):
       """returns -self."""
       return self*(-1.)

   def __pos__(self):
       """returns +self."""
       return self

   def __abs__(self):
       """returns absolute value"""
       return Abs_Symbol(self)

   def __sub__(self,other):
       """subtracts other from symbol self. if _testForZero(other) self is returned."""
       if _testForZero(other):
          return self
       else:
          return self+(-other)

   def __rsub__(self,other):
       """subtracts symbol self from other."""
       return -self+other

   def __div__(self,other):
       """divides symbol self by other."""
       if isinstance(other,Symbol):
          a=_matchShape([self,other])
          return Div_Symbol(a[0],a[1])
       else:
          return self*(1./other)

   def __rdiv__(self,other):
       """dived other by symbol self. if _testForZero(other) 0 is returned."""
       if _testForZero(other):
          return 0
       else:
          a=_matchShape([self,other])
          return Div_Symbol(a[0],a[1])

   def __pow__(self,other):
       """raises symbol self to the power of other"""
       a=_matchShape([self,other])
       return Power_Symbol(a[0],a[1])

   def __rpow__(self,other):
       """raises other to the symbol self"""
       a=_matchShape([self,other])
       return Power_Symbol(a[1],a[0])

   def __mul__(self,other):
       """multiplies other by symbol self. if _testForZero(other) 0 is returned."""
       if _testForZero(other):
          return 0
       else:
          a=_matchShape([self,other])
          return Mult_Symbol(a[0],a[1])

   def __rmul__(self,other):
       """multiplies other by symbol self. if _testSForZero(other) 0 is returned."""
       return self*other

   def __getitem__(self,sl):
          print sl

def Float_Symbol(Symbol):
    def __init__(self,name="symbol",shape=(),args=[]):
        Symbol.__init__(self,dim=0,name="symbol",shape=(),args=[])

class ScalarSymbol(Symbol):
   """a scalar symbol"""
   def __init__(self,dim=3,name="scalar"):
      """creates a scalar symbol of spatial dimension dim"""
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(),dim=d,name=name)

class VectorSymbol(Symbol):
   """a vector symbol"""
   def __init__(self,dim=3,name="vector"):
      """creates a vector symbol of spatial dimension dim"""
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,),dim=d,name=name)

class TensorSymbol(Symbol):
   """a tensor symbol""" 
   def __init__(self,dim=3,name="tensor"):
      """creates a tensor symbol of spatial dimension dim"""
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,d),dim=d,name=name)

class Tensor3Symbol(Symbol):
   """a tensor order 3 symbol"""
   def __init__(self,dim=3,name="tensor3"):
      """creates a tensor order 3 symbol of spatial dimension dim"""
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,d,d),dim=d,name=name)

class Tensor4Symbol(Symbol):
   """a tensor order 4 symbol"""
   def __init__(self,dim=3,name="tensor4"):
      """creates a tensor order 4 symbol of spatial dimension dim"""    
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,d,d,d),dim=d,name=name)

class Add_Symbol(Symbol):
   """symbol representing the sum of two arguments"""
   def __init__(self,arg0,arg1):
       a=[arg0,arg1]
       Symbol.__init__(self,dim=_extractDim(a),shape=_extractShape(a),args=a)
   def __str__(self):
      return "(%s+%s)"%(str(self.getArgument(0)),str(self.getArgument(1)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return a[0]+a[1]
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return a[0]+a[1]

class Mult_Symbol(Symbol):
   """symbol representing the product of two arguments"""
   def __init__(self,arg0,arg1):
       a=[arg0,arg1]
       Symbol.__init__(self,dim=_extractDim(a),shape=_extractShape(a),args=a)
   def __str__(self):
      return "(%s*%s)"%(str(self.getArgument(0)),str(self.getArgument(1)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return a[0]*a[1]
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return self.getArgument(1)*a[0]+self.getArgument(0)*a[1]

class Div_Symbol(Symbol):
   """symbol representing the quotient of two arguments"""
   def __init__(self,arg0,arg1):
       a=[arg0,arg1]
       Symbol.__init__(self,dim=_extractDim(a),shape=_extractShape(a),args=a)
   def __str__(self):
      return "(%s/%s)"%(str(self.getArgument(0)),str(self.getArgument(1)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return a[0]/a[1]
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return (a[0]*self.getArgument(1)-self.getArgument(0)*a[1])/ \
                          (self.getArgument(1)*self.getArgument(1))

class Power_Symbol(Symbol):
   """symbol representing the power of the first argument to the power of the second argument"""
   def __init__(self,arg0,arg1):
       a=[arg0,arg1]
       Symbol.__init__(self,dim=_extractDim(a),shape=_extractShape(a),args=a)
   def __str__(self):
      return "(%s**%s)"%(str(self.getArgument(0)),str(self.getArgument(1)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return a[0]**a[1]
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return self*(a[1]*log(self.getArgument(0))+self.getArgument(1)/self.getArgument(0)*a[0])

class Abs_Symbol(Symbol):
   """symbol representing absolute value of its argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "abs(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return abs(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return sign(self.getArgument(0))*self.getDifferentiatedArguments(arg)[0]

#=========================================================
#   some little helpers
#=========================================================
def _testForZero(arg):
   """returns True is arg is considered of being zero"""
   if isinstance(arg,int):
      return not arg>0
   elif isinstance(arg,float):
      return not arg>0.
   elif isinstance(arg,numarray.NumArray): 
      a=abs(arg)
      while isinstance(a,numarray.NumArray): a=numarray.sometrue(a)
      return not a>0
   else:
      return False

def _extractDim(args):
    dim=None
    for a in args:
       if hasattr(a,"getDim"):
          d=a.getDim()
          if dim==None: 
             dim=d
          else:
             if dim!=d: raise ValueError,"inconsistent spatial dimension of arguments"
    if dim==None:
       raise ValueError,"cannot recover spatial dimension"
    return dim

def _identifyShape(arg):
   """identifies the shape of arg."""
   if hasattr(arg,"getShape"):
       arg_shape=arg.getShape()
   elif hasattr(arg,"shape"):
     s=arg.shape
     if callable(s):
       arg_shape=s()
     else:
       arg_shape=s
   else:
       arg_shape=()
   return arg_shape

def _extractShape(args):
    """extracts the common shape of the list of arguments args"""
    shape=None
    for a in args:
       a_shape=_identifyShape(a)
       if shape==None: shape=a_shape
       if shape!=a_shape: raise ValueError,"inconsistent shape"
    if shape==None:
       raise ValueError,"cannot recover shape"
    return shape

def _matchShape(args,shape=None):
    """returns the list of arguments args as object which have all the specified shape.
       if shape is not given the shape "largest" shape of args is used."""
    # identify the list of shapes:
    arg_shapes=[]
    for a in args: arg_shapes.append(_identifyShape(a))
    # get the largest shape (currently the longest shape):
    if shape==None: shape=max(arg_shapes)
    
    out=[]
    for i in range(len(args)):
       if shape==arg_shapes[i]:
          out.append(args[i])
       else:
          if len(shape)==0: # then len(arg_shapes[i])>0
            raise ValueError,"cannot adopt shape of %s to %s"%(str(args[i]),str(shape))
          else: 
            if len(arg_shapes[i])==0:
                out.append(outer(args[i],numarray.ones(shape)))         
            else:  
                raise ValueError,"cannot adopt shape of %s to %s"%(str(args[i]),str(shape))
    return out  

#=========================================================
#   wrappers for various mathematical functions:
#=========================================================
def diff(arg,dep):
    """returns the derivative of arg with respect to dep. If arg is not Symbol object
       0 is returned"""
    if isinstance(arg,Symbol):
       return arg.diff(dep)
    elif hasattr(arg,"shape"):
          if callable(arg.shape):
              return numarray.zeros(arg.shape(),numarray.Float)
          else:
              return numarray.zeros(arg.shape,numarray.Float)
    else: 
       return 0

def exp(arg):
    """
    @brief applies the exponential function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Exp_Symbol(arg)
    elif hasattr(arg,"exp"):
       return arg.exp()
    else:
       return numarray.exp(arg)

class Exp_Symbol(Symbol):
   """symbol representing the power of the first argument to the power of the second argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "exp(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return exp(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return self*self.getDifferentiatedArguments(arg)[0]

def sqrt(arg):
    """
    @brief applies the squre root function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Sqrt_Symbol(arg)
    elif hasattr(arg,"sqrt"):
       return arg.sqrt()
    else:
       return numarray.sqrt(arg)       

class Sqrt_Symbol(Symbol):
   """symbol representing square root of argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "sqrt(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return sqrt(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return (-0.5)/self*self.getDifferentiatedArguments(arg)[0]

def log(arg):
    """
    @brief applies the logarithmic function bases exp(1.) to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Log_Symbol(arg)
    elif hasattr(arg,"log"):
       return arg.log()
    else:
       return numarray.log(arg)

class Log_Symbol(Symbol):
   """symbol representing logarithm of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "log(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return log(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return self.getDifferentiatedArguments(arg)[0]/self.getArgument(0)

def ln(arg):
    """
    @brief applies the natural logarithmic function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Ln_Symbol(arg)
    elif hasattr(arg,"ln"):
       return arg.log()
    else:
       return numarray.log(arg)

class Ln_Symbol(Symbol):
   """symbol representing natural logarithm of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "ln(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return ln(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return self.getDifferentiatedArguments(arg)[0]/self.getArgument(0)

def sin(arg):
    """
    @brief applies the sinus function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Sin_Symbol(arg)
    elif hasattr(arg,"sin"):
       return arg.sin()
    else:
       return numarray.sin(arg)

class Sin_Symbol(Symbol):
   """symbol representing logarithm of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "sin(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return sin(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return cos(self.getArgument(0))*self.getDifferentiatedArguments(arg)[0]

def cos(arg):
    """
    @brief applies the sinus function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Cos_Symbol(arg)
    elif hasattr(arg,"cos"):
       return arg.cos()
    else:
       return numarray.cos(arg)

class Cos_Symbol(Symbol):
   """symbol representing logarithm of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "cos(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return cos(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return -sin(self.getArgument(0))*self.getDifferentiatedArguments(arg)[0]

def tan(arg):
    """
    @brief applies the sinus function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Tan_Symbol(arg)
    elif hasattr(arg,"tan"):
       return arg.tan()
    else:
       return numarray.tan(arg)

class Tan_Symbol(Symbol):
   """symbol representing logarithm of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "tan(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return tan(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       s=cos(self.getArgument(0))
       return 1./(s*s)*self.getDifferentiatedArguments(arg)[0]

def sign(arg):
    """
    @brief applies the sign function to arg
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Sign_Symbol(arg)
    elif hasattr(arg,"sign"):
       return arg.sign()
    else:
       return numarray.greater(arg,numarray.zeros(arg.shape,numarray.Float))- \
              numarray.less(arg,numarray.zeros(arg.shape,numarray.Float))

class Sign_Symbol(Symbol):
   """symbol representing the sign of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "sign(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return sign(self.getEvaluatedArguments(argval)[0])

def maxval(arg):
    """
    @brief returns the maximum value of argument arg""
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Max_Symbol(arg)
    elif hasattr(arg,"maxval"):
       return arg.maxval()
    elif hasattr(arg,"max"):
       return arg.max()
    else:
       return arg

class Max_Symbol(Symbol):
   """symbol representing the sign of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "maxval(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return maxval(self.getEvaluatedArguments(argval)[0])

def minval(arg):
    """
    @brief returns the maximum value of argument arg""
    @param arg (input): argument 
    """
    if isinstance(arg,Symbol):
       return Min_Symbol(arg)
    elif hasattr(arg,"maxval"):
       return arg.minval()
    elif hasattr(arg,"min"):
       return arg.min()
    else:
       return arg

class Min_Symbol(Symbol):
   """symbol representing the sign of the argument"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "minval(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return minval(self.getEvaluatedArguments(argval)[0])

def wherePositive(arg):
    """
    @brief returns the maximum value of argument arg""
    @param arg (input): argument 
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,Symbol):
       return WherePositive_Symbol(arg)
    elif hasattr(arg,"wherePositive"):
       return arg.minval()
    elif hasattr(arg,"wherePositive"):
       numarray.greater(arg,numarray.zeros(arg.shape,numarray.Float))
    else:
       if arg>0:
          return 1.
       else:
          return 0.

class WherePositive_Symbol(Symbol):
   """symbol representing the wherePositive function"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "wherePositive(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return wherePositive(self.getEvaluatedArguments(argval)[0])

def whereNegative(arg):
    """
    @brief returns the maximum value of argument arg""
    @param arg (input): argument 
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,Symbol):
       return WhereNegative_Symbol(arg)
    elif hasattr(arg,"whereNegative"):
       return arg.whereNegative()
    elif hasattr(arg,"shape"):
       numarray.less(arg,numarray.zeros(arg.shape,numarray.Float))
    else:
       if arg<0:
          return 1.
       else:
          return 0.

class WhereNegative_Symbol(Symbol):
   """symbol representing the whereNegative function"""
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "whereNegative(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return whereNegative(self.getEvaluatedArguments(argval)[0])

def outer(arg0,arg1):
   if _testForZero(arg0) or _testForZero(arg1):
      return 0
   else:
      if isinstance(arg0,Symbol) or isinstance(arg1,Symbol):
        return Outer_Symbol(arg0,arg1)
      elif _identifyShape(arg0)==() or _identifyShape(arg1)==():
        return arg0*arg1
      elif isinstance(arg0,numarray.NumArray) and isinstance(arg1,numarray.NumArray): 
        return numarray.outer(arg0,arg1)
      else:
        if arg0.getRank()==1 and arg1.getRank()==1:
          out=escript.Data(0,(arg0.getShape()[0],arg1.getShape()[0]),arg1.getFunctionSpace())
          for i in range(arg0.getShape()[0]):
            for j in range(arg1.getShape()[0]):
                out[i,j]=arg0[i]*arg1[j]
          return out
        else:
          raise ValueError,"outer is not fully implemented yet."

class Outer_Symbol(Symbol):
   """symbol representing the outer product of its two argument"""
   def __init__(self,arg0,arg1):
       a=[arg0,arg1]
       s=tuple(list(_identifyShape(arg0))+list(_identifyShape(arg1)))
       Symbol.__init__(self,shape=s,dim=_extractDim(a),args=a)
   def __str__(self):
      return "outer(%s,%s)"%(str(self.getArgument(0)),str(self.getArgument(1)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return outer(a[0],a[1])
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return outer(a[0],self.getArgument(1))+outer(self.getArgument(0),a[1])

def interpolate(arg,where):
    """
    @brief interpolates the function into the FunctionSpace where.

    @param arg    interpolant 
    @param where  FunctionSpace to interpolate to
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,Symbol):
       return Interpolated_Symbol(arg,where)
    else:
       return escript.Data(arg,where)

def Interpolated_Symbol(Symbol):
   """symbol representing the integral of the argument"""
   def __init__(self,arg,where):
        Symbol.__init__(self,shape=_extractShape(arg),dim=_extractDim([arg]),args=[arg,where])
   def __str__(self):
      return "interpolated(%s)"%(str(self.getArgument(0)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return integrate(a[0],where=self.getArgument(1))
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return integrate(a[0],where=self.getArgument(1))

def grad(arg,where=None):
    """
    @brief returns the spatial gradient of arg at where.

    @param arg:   Data object representing the function which gradient to be calculated.
    @param where: FunctionSpace in which the gradient will be calculated. If not present or
                  None an appropriate default is used. 
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,Symbol):
       return Grad_Symbol(arg,where)
    elif hasattr(arg,"grad"):
       if where==None:
          return arg.grad()
       else:
          return arg.grad(where)
    else:
       return arg*0.

def Grad_Symbol(Symbol):
   """symbol representing the gradient of the argument"""
   def __init__(self,arg,where=None):
       d=_extractDim([arg])
       s=tuple(list(_identifyShape([arg])).append(d))
       Symbol.__init__(self,shape=s,dim=_extractDim([arg]),args=[arg,where])
   def __str__(self):
      return "grad(%s)"%(str(self.getArgument(0)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return grad(a[0],where=self.getArgument(1))
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return grad(a[0],where=self.getArgument(1))

def integrate(arg,where=None):
    """
    @brief return the integral if the function represented by Data object arg over its domain.

    @param arg:   Data object representing the function which is integrated.
    @param where: FunctionSpace in which the integral is calculated. If not present or
                  None an appropriate default is used.
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,Symbol):
       return Integral_Symbol(arg,where)
    else:    
       if not where==None: arg=escript.Data(arg,where)
       if arg.getRank()==0:
         return arg.integrate()[0]
       else:
         return arg.integrate()

def Integral_Symbol(Float_Symbol):
   """symbol representing the integral of the argument"""
   def __init__(self,arg,where=None):
       Float_Symbol.__init__(self,shape=_identifyShape([arg]),args=[arg,where])
   def __str__(self):
      return "integral(%s)"%(str(self.getArgument(0)))
   def eval(self,argval):
       a=self.getEvaluatedArguments(argval)
       return integrate(a[0],where=self.getArgument(1))
   def _diff(self,arg):
       a=self.getDifferentiatedArguments(arg)
       return integrate(a[0],where=self.getArgument(1))

#=============================
#
# wrapper for various functions: if the argument has attribute the function name
# as an argument it calls the corresponding methods. Otherwise the corresponding
# numarray function is called.

# functions involving the underlying Domain:


# functions returning Data objects:

def transpose(arg,axis=None):
    """
    @brief returns the transpose of the Data object arg. 

    @param arg
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
    @brief return 

    @param arg
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
       for i in range(s[0]):
          for j in range(s[1]):
             out+=arg[i,j]
       return out
       # end hack for transpose 
       return arg.transpose(axis0=axis0,axis1=axis1)
    else:
       return numarray.trace(arg,axis0=axis0,axis1=axis1)

def Trace_Symbol(Symbol):
    pass

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

    @param arg0
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

def inner(arg0,arg1):
    """
    @brief

    @param arg0, arg1
    """
    sum=escript.Scalar(0,arg0.getFunctionSpace())
    if arg.getRank()==0:
          return arg0*arg1
    elif arg.getRank()==1:
         sum=escript.Scalar(0,arg.getFunctionSpace())
         for i in range(arg.getShape()[0]):
            sum+=arg0[i]*arg1[i]
    elif arg.getRank()==2:
        sum=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
           for j in range(arg.getShape()[1]):
              sum+=arg0[i,j]*arg1[i,j]
    elif arg.getRank()==3:
        sum=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
            for j in range(arg.getShape()[1]):
               for k in range(arg.getShape()[2]):
                  sum+=arg0[i,j,k]*arg1[i,j,k]
    elif arg.getRank()==4:
        sum=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
           for j in range(arg.getShape()[1]):
              for k in range(arg.getShape()[2]):
                 for l in range(arg.getShape()[3]):
                    sum+=arg0[i,j,k,l]*arg1[i,j,k,l]
    else:
          raise SystemError,"inner is not been implemented yet"
    return sum

def matrixmult(arg0,arg1):

    if isinstance(arg1,numarray.NumArray) and isinstance(arg0,numarray.NumArray):
        numarray.matrixmult(arg0,arg1)
    else:
      # escript.matmult(arg0,arg1)
      if isinstance(arg1,escript.Data) and not isinstance(arg0,escript.Data):
        arg0=escript.Data(arg0,arg1.getFunctionSpace())
      elif isinstance(arg0,escript.Data) and not isinstance(arg1,escript.Data):
        arg1=escript.Data(arg1,arg0.getFunctionSpace())
      if arg0.getRank()==2 and arg1.getRank()==1:
          out=escript.Data(0,(arg0.getShape()[0],),arg0.getFunctionSpace())
          for i in range(arg0.getShape()[0]):
             for j in range(arg0.getShape()[1]):
               out[i]+=arg0[i,j]*arg1[j]
          return out
      else:
          raise SystemError,"matrixmult is not fully implemented yet!"

#=========================================================
# reduction operations:
#=========================================================
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

def dot(arg0,arg1):
    """
    @brief

    @param arg
    """
    if isinstance(arg0,escript.Data):
       return arg0.dot(arg1)
    elif isinstance(arg1,escript.Data):
       return arg1.dot(arg0)
    else:
       return numarray.dot(arg0,arg1)

def kronecker(d):
   if hasattr(d,"getDim"):
      return numarray.identity(d.getDim())
   else:
      return numarray.identity(d)

def unit(i,d):
   """
   @brief return a unit vector of dimension d with nonzero index i
   @param d dimension
   @param i index
   """
   e = numarray.zeros((d,),numarray.Float)
   e[i] = 1.0
   return e

# ============================================
#   testing
# ============================================

if __name__=="__main__":
  u=ScalarSymbol(dim=2,name="u")
  v=ScalarSymbol(dim=2,name="v")
  v=VectorSymbol(2,"v")
  u=VectorSymbol(2,"u")

  print u+5,(u+5).diff(u)
  print 5+u,(5+u).diff(u)
  print u+v,(u+v).diff(u)
  print v+u,(v+u).diff(u)

  print u*5,(u*5).diff(u)
  print 5*u,(5*u).diff(u)
  print u*v,(u*v).diff(u)
  print v*u,(v*u).diff(u)

  print u-5,(u-5).diff(u)
  print 5-u,(5-u).diff(u) 
  print u-v,(u-v).diff(u)
  print v-u,(v-u).diff(u)

  print u/5,(u/5).diff(u)
  print 5/u,(5/u).diff(u)
  print u/v,(u/v).diff(u)
  print v/u,(v/u).diff(u)

  print u**5,(u**5).diff(u)
  print 5**u,(5**u).diff(u)
  print u**v,(u**v).diff(u)
  print v**u,(v**u).diff(u)

  print exp(u),exp(u).diff(u)
  print sqrt(u),sqrt(u).diff(u)
  print log(u),log(u).diff(u)
  print sin(u),sin(u).diff(u)
  print cos(u),cos(u).diff(u)
  print tan(u),tan(u).diff(u)
  print sign(u),sign(u).diff(u)
  print abs(u),abs(u).diff(u)
  print wherePositive(u),wherePositive(u).diff(u)
  print whereNegative(u),whereNegative(u).diff(u)
  print maxval(u),maxval(u).diff(u)
  print minval(u),minval(u).diff(u)

  g=grad(u)
  print diff(5*g,g)
  4*(g+transpose(g))/2+6*trace(g)*kronecker(3) 

#
# $Log$
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
