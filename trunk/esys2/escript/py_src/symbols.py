# $Id$

"""
symbolic tool box for escript
"""

import numarray

#===========================================================
# a simple tool box to deal with _differentials of functions 
#===========================================================

class Symbol:
   """
   Symbol class.
   """
   def __init__(self,name="symbol",shape=(),dim=3,args=[]):
       """
       Creates an instance of a symbol of shape shape and spatial dimension
       dim.
       
       The symbol may depending on a list of arguments args which may be
       symbols or other objects. name gives the name of the symbol.
       """

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
       """
       Returns the i-th argument.
       """
       return self.__args[i]

   def getDim(self):
       """
       Returns the spatial dimension of the symbol.
       """
       return self.__dim

   def getRank(self):
       """
       Returns the rank of the symbol.
       """
       return len(self.getShape())

   def getShape(self):
       """
       Returns the shape of the symbol.
       """
       return self.__shape

   def getEvaluatedArguments(self,argval):
       """
       Returns the list of evaluated arguments by subsituting symbol u by 
       argval[u].
       """
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
       """
       Returns the list of the arguments _differentiated by arg.
       """ 
       out=[]
       for a in self.__args:
          if isinstance(a,Symbol):
            out.append(a.diff(arg))
          else:
            out.append(0)
       return out

   def diff(self,arg):
       """
       Returns the _differention of self by arg.
       """
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
       """
       Return derivate of self with respect to arg (!=self). 

       This method is overwritten by a particular symbol.
       """
       return 0

   def eval(self,argval):
       """
       Subsitutes symbol u in self by argval[u] and returns the result. If
       self is not a key of argval then self is returned.
       """
       if argval.has_key(self):
         return argval[self]
       else:
         return self

   def __str__(self):
       """
       Returns a string representation of the symbol.
       """
       return self.__name

   def __add__(self,other):
       """
       Adds other to symbol self. if _testForZero(other) self is returned.
       """
       if _testForZero(other):
          return self
       else:
          a=_matchShape([self,other])
          return Add_Symbol(a[0],a[1])

   def __radd__(self,other):
       """
       Adds other to symbol self. if _testForZero(other) self is returned.
       """
       return self+other

   def __neg__(self):
       """
       Returns -self.
       """
       return self*(-1.)

   def __pos__(self):
       """
       Returns +self.
       """
       return self

   def __abs__(self):
       """
       Returns absolute value.
       """
       return Abs_Symbol(self)

   def __sub__(self,other):
       """
       Subtracts other from symbol self. 
       
       If _testForZero(other) self is returned.
       """
       if _testForZero(other):
          return self
       else:
          return self+(-other)

   def __rsub__(self,other):
       """
       Subtracts symbol self from other.
       """
       return -self+other

   def __div__(self,other):
       """
       Divides symbol self by other.
       """
       if isinstance(other,Symbol):
          a=_matchShape([self,other])
          return Div_Symbol(a[0],a[1])
       else:
          return self*(1./other)

   def __rdiv__(self,other):
       """
       Dived other by symbol self. if _testForZero(other) 0 is returned.
       """
       if _testForZero(other):
          return 0
       else:
          a=_matchShape([self,other])
          return Div_Symbol(a[0],a[1])

   def __pow__(self,other):
       """
       Raises symbol self to the power of other.
       """
       a=_matchShape([self,other])
       return Power_Symbol(a[0],a[1])

   def __rpow__(self,other):
       """
       Raises other to the symbol self.
       """
       a=_matchShape([self,other])
       return Power_Symbol(a[1],a[0])

   def __mul__(self,other):
       """
       Multiplies other by symbol self. if _testForZero(other) 0 is returned.
       """
       if _testForZero(other):
          return 0
       else:
          a=_matchShape([self,other])
          return Mult_Symbol(a[0],a[1])

   def __rmul__(self,other):
       """
       Multiplies other by symbol self. if _testSForZero(other) 0 is returned.
       """
       return self*other

   def __getitem__(self,sl):
          print sl

class Float_Symbol(Symbol):
    def __init__(self,name="symbol",shape=(),args=[]):
        Symbol.__init__(self,dim=0,name="symbol",shape=(),args=[])

class ScalarSymbol(Symbol):
   """
   A scalar symbol.
   """
   def __init__(self,dim=3,name="scalar"):
      """
      Creates a scalar symbol of spatial dimension dim.
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(),dim=d,name=name)

class VectorSymbol(Symbol):
   """
   A vector symbol.
   """
   def __init__(self,dim=3,name="vector"):
      """
      Creates a vector symbol of spatial dimension dim.
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,),dim=d,name=name)

class TensorSymbol(Symbol):
   """
   A tensor symbol.
   """ 
   def __init__(self,dim=3,name="tensor"):
      """
      Creates a tensor symbol of spatial dimension dim.
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,d),dim=d,name=name)

class Tensor3Symbol(Symbol):
   """
   A tensor order 3 symbol.
   """
   def __init__(self,dim=3,name="tensor3"):
      """
      Creates a tensor order 3 symbol of spatial dimension dim.
      """
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,d,d),dim=d,name=name)

class Tensor4Symbol(Symbol):
   """
   A tensor order 4 symbol.
   """
   def __init__(self,dim=3,name="tensor4"):
      """
      Creates a tensor order 4 symbol of spatial dimension dim.
      """ 
      if hasattr(dim,"getDim"):
           d=dim.getDim()
      else:    
           d=dim
      Symbol.__init__(self,shape=(d,d,d,d),dim=d,name=name)

class Add_Symbol(Symbol):
   """
   Symbol representing the sum of two arguments.
   """
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
   """
   Symbol representing the product of two arguments.
   """
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
   """
   Symbol representing the quotient of two arguments.
   """
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
   """
   Symbol representing the power of the first argument to the power of the
   second argument.
   """
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
   """
   Symbol representing absolute value of its argument.
   """
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
   """
   Returns True is arg is considered to be zero.
   """
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
   """
   Identifies the shape of arg.
   """
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
    """
    Extracts the common shape of the list of arguments args.
    """
    shape=None
    for a in args:
       a_shape=_identifyShape(a)
       if shape==None: shape=a_shape
       if shape!=a_shape: raise ValueError,"inconsistent shape"
    if shape==None:
       raise ValueError,"cannot recover shape"
    return shape

def _matchShape(args,shape=None):
    """
    Returns the list of arguments args as object which have all the 
    specified shape.

    If shape is not given the shape "largest" shape of args is used.
    """
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

class Exp_Symbol(Symbol):
   """
   Symbol representing the power of the first argument to the power of the
   second argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "exp(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return exp(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return self*self.getDifferentiatedArguments(arg)[0]

class Sqrt_Symbol(Symbol):
   """
   Symbol representing square root of argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "sqrt(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return sqrt(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return (-0.5)/self*self.getDifferentiatedArguments(arg)[0]

class Log_Symbol(Symbol):
   """
   Symbol representing logarithm of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "log(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return log(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return self.getDifferentiatedArguments(arg)[0]/self.getArgument(0)

class Ln_Symbol(Symbol):
   """
   Symbol representing natural logarithm of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "ln(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return ln(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return self.getDifferentiatedArguments(arg)[0]/self.getArgument(0)

class Sin_Symbol(Symbol):
   """
   Symbol representing sin of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "sin(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return sin(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return cos(self.getArgument(0))*self.getDifferentiatedArguments(arg)[0]

class Cos_Symbol(Symbol):
   """
   Symbol representing cos of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "cos(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return cos(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       return -sin(self.getArgument(0))*self.getDifferentiatedArguments(arg)[0]

class Tan_Symbol(Symbol):
   """
   Symbol representing tan of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "tan(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return tan(self.getEvaluatedArguments(argval)[0])
   def _diff(self,arg):
       s=cos(self.getArgument(0))
       return 1./(s*s)*self.getDifferentiatedArguments(arg)[0]

class Sign_Symbol(Symbol):
   """
   Symbol representing the sign of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "sign(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return sign(self.getEvaluatedArguments(argval)[0])

class Max_Symbol(Symbol):
   """
   Symbol representing the maximum value of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "maxval(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return maxval(self.getEvaluatedArguments(argval)[0])

class Min_Symbol(Symbol):
   """
   Symbol representing the minimum value of the argument.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "minval(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return minval(self.getEvaluatedArguments(argval)[0])

class WherePositive_Symbol(Symbol):
   """
   Symbol representing the wherePositive function.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "wherePositive(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return wherePositive(self.getEvaluatedArguments(argval)[0])

class WhereNegative_Symbol(Symbol):
   """
   Symbol representing the whereNegative function.
   """
   def __init__(self,arg):
       Symbol.__init__(self,shape=arg.getShape(),dim=arg.getDim(),args=[arg])
   def __str__(self):
      return "whereNegative(%s)"%str(self.getArgument(0))
   def eval(self,argval):
       return whereNegative(self.getEvaluatedArguments(argval)[0])

class Outer_Symbol(Symbol):
   """
   Symbol representing the outer product of its two arguments.
   """
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

class Interpolated_Symbol(Symbol):
   """
   Symbol representing the integral of the argument.
   """
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

class Grad_Symbol(Symbol):
   """
   Symbol representing the gradient of the argument.
   """
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

class Integral_Symbol(Float_Symbol):
   """
   Symbol representing the integral of the argument.
   """
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
# Revision 1.2  2005/09/15 03:44:19  jgs
# Merge of development branch dev-02 back to main trunk on 2005-09-15
#
# Revision 1.1.2.1  2005/09/07 10:32:05  gross
# Symbols removed from util and put into symmbols.py.
#
#
