# $Id$

## @file util.py

"""
Utility functions for escript

@todo:

  - binary operations @ (@=+,-,*,/,**)::
      (a@b)[:,*]=a[:]@b[:,*] if rank(a)<rank(b)
      (a@b)[:]=a[:]@b[:] if rank(a)=rank(b)
      (a@b)[*,:]=a[*,:]@b[:] if rank(a)>rank(b)
  - implementation of outer::
      outer(a,b)[:,*]=a[:]*b[*]
  - trace:: 
      trace(arg,axis0=a0,axis1=a1)(:,&,*)=sum_i trace(:,i,&,i,*) (i are at index a0 and a1)
"""

import numarray
import escript
import symbols
import os

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

#=========================================================
def saveVTK(filename,domain=None,**data):
    """
    writes a L{Data} objects into a files using the the VTK XML file format.

    Example:

       tmp=Scalar(..)
       v=Vector(..)
       saveVTK("solution.xml",temperature=tmp,velovity=v)

    tmp and v are written into "solution.dx" where tmp is named "temperature" and v is named "velovity"

    @param filename: file name of the output file
    @type filename: C(str}
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

#=========================================================
def saveDX(filename,domain=None,**data):
    """
    writes a L{Data} objects into a files using the the DX file format.

    Example:

       tmp=Scalar(..)
       v=Vector(..)
       saveDX("solution.dx",temperature=tmp,velovity=v)

    tmp and v are written into "solution.dx" where tmp is named "temperature" and v is named "velovity". 

    @param filename: file name of the output file
    @type filename: C(str}
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

#=========================================================

def exp(arg):
    """
    Applies the exponential function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Exp_Symbol(arg)
    elif hasattr(arg,"exp"):
       return arg.exp()
    else:
       return numarray.exp(arg)

def sqrt(arg):
    """
    Applies the squre root function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Sqrt_Symbol(arg)
    elif hasattr(arg,"sqrt"):
       return arg.sqrt()
    else:
       return numarray.sqrt(arg)       

def log(arg):
    """
    Applies the logarithmic function base 10 to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Log_Symbol(arg)
    elif hasattr(arg,"log"):
       return arg.log()
    else:
       return numarray.log10(arg)

def ln(arg):
    """
    Applies the natural logarithmic function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Ln_Symbol(arg)
    elif hasattr(arg,"ln"):
       return arg.ln()
    else:
       return numarray.log(arg)

def sin(arg):
    """
    Applies the sin function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Sin_Symbol(arg)
    elif hasattr(arg,"sin"):
       return arg.sin()
    else:
       return numarray.sin(arg)

def cos(arg):
    """
    Applies the cos function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Cos_Symbol(arg)
    elif hasattr(arg,"cos"):
       return arg.cos()
    else:
       return numarray.cos(arg)

def tan(arg):
    """
    Applies the tan function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Tan_Symbol(arg)
    elif hasattr(arg,"tan"):
       return arg.tan()
    else:
       return numarray.tan(arg)

def asin(arg):
    """
    Applies the asin function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Asin_Symbol(arg)
    elif hasattr(arg,"asin"):
       return arg.asin()
    else:
       return numarray.asin(arg)

def acos(arg):
    """
    Applies the acos function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Acos_Symbol(arg)
    elif hasattr(arg,"acos"):
       return arg.acos()
    else:
       return numarray.acos(arg)

def atan(arg):
    """
    Applies the atan function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Atan_Symbol(arg)
    elif hasattr(arg,"atan"):
       return arg.atan()
    else:
       return numarray.atan(arg)

def sinh(arg):
    """
    Applies the sinh function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Sinh_Symbol(arg)
    elif hasattr(arg,"sinh"):
       return arg.sinh()
    else:
       return numarray.sinh(arg)

def cosh(arg):
    """
    Applies the cosh function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Cosh_Symbol(arg)
    elif hasattr(arg,"cosh"):
       return arg.cosh()
    else:
       return numarray.cosh(arg)

def tanh(arg):
    """
    Applies the tanh function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Tanh_Symbol(arg)
    elif hasattr(arg,"tanh"):
       return arg.tanh()
    else:
       return numarray.tanh(arg)

def asinh(arg):
    """
    Applies the asinh function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Asinh_Symbol(arg)
    elif hasattr(arg,"asinh"):
       return arg.asinh()
    else:
       return numarray.asinh(arg)

def acosh(arg):
    """
    Applies the acosh function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Acosh_Symbol(arg)
    elif hasattr(arg,"acosh"):
       return arg.acosh()
    else:
       return numarray.acosh(arg)

def atanh(arg):
    """
    Applies the atanh function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Atanh_Symbol(arg)
    elif hasattr(arg,"atanh"):
       return arg.atanh()
    else:
       return numarray.atanh(arg)

def sign(arg):
    """
    Applies the sign function to arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Sign_Symbol(arg)
    elif hasattr(arg,"sign"):
       return arg.sign()
    else:
       return numarray.greater(arg,numarray.zeros(arg.shape,numarray.Float))- \
              numarray.less(arg,numarray.zeros(arg.shape,numarray.Float))

def maxval(arg):
    """
    Returns the maximum value of argument arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Max_Symbol(arg)
    elif hasattr(arg,"maxval"):
       return arg.maxval()
    elif hasattr(arg,"max"):
       return arg.max()
    else:
       return arg

def minval(arg):
    """
    Returns the minimum value of argument arg.

    @param arg: argument 
    """
    if isinstance(arg,symbols.Symbol):
       return symbols.Min_Symbol(arg)
    elif hasattr(arg,"maxval"):
       return arg.minval()
    elif hasattr(arg,"min"):
       return arg.min()
    else:
       return arg

def wherePositive(arg):
    """
    Returns the positive values of argument arg.

    @param arg: argument 
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,symbols.Symbol):
       return symbols.WherePositive_Symbol(arg)
    elif hasattr(arg,"wherePositive"):
       return arg.minval()
    elif hasattr(arg,"wherePositive"):
       numarray.greater(arg,numarray.zeros(arg.shape,numarray.Float))
    else:
       if arg>0:
          return 1.
       else:
          return 0.

def whereNegative(arg):
    """
    Returns the negative values of argument arg.

    @param arg: argument 
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,symbols.Symbol):
       return symbols.WhereNegative_Symbol(arg)
    elif hasattr(arg,"whereNegative"):
       return arg.whereNegative()
    elif hasattr(arg,"shape"):
       numarray.less(arg,numarray.zeros(arg.shape,numarray.Float))
    else:
       if arg<0:
          return 1.
       else:
          return 0.

def maximum(arg0,arg1):
    """
    Return arg1 where arg1 is bigger then arg0 otherwise arg0 is returned.
    """
    m=whereNegative(arg0-arg1)
    return m*arg1+(1.-m)*arg0
   
def minimum(arg0,arg1):
    """
    Return arg0 where arg1 is bigger then arg0 otherwise arg1 is returned.
    """
    m=whereNegative(arg0-arg1)
    return m*arg0+(1.-m)*arg1
   
def outer(arg0,arg1):
   if _testForZero(arg0) or _testForZero(arg1):
      return 0
   else:
      if isinstance(arg0,symbols.Symbol) or isinstance(arg1,symbols.Symbol):
        return symbols.Outer_Symbol(arg0,arg1)
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

def interpolate(arg,where):
    """
    Interpolates the function into the FunctionSpace where.

    @param arg:    interpolant 
    @param where:  FunctionSpace to interpolate to
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,symbols.Symbol):
       return symbols.Interpolated_Symbol(arg,where)
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
    return trace(grad(arg,where))

def jump(arg):
    """
    Returns the jump of arg across a continuity.

    @param arg:   Data object representing the function which gradient 
                  to be calculated.
    """
    d=arg.getDomain()
    return arg.interpolate(escript.FunctionOnContactOne())-arg.interpolate(escript.FunctionOnContactZero())
   

def grad(arg,where=None):
    """
    Returns the spatial gradient of arg at where.

    @param arg:   Data object representing the function which gradient 
                  to be calculated.
    @param where: FunctionSpace in which the gradient will be calculated. 
                  If not present or C{None} an appropriate default is used. 
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,symbols.Symbol):
       return symbols.Grad_Symbol(arg,where)
    elif hasattr(arg,"grad"):
       if where==None:
          return arg.grad()
       else:
          return arg.grad(where)
    else:
       return arg*0.

def integrate(arg,where=None):
    """
    Return the integral if the function represented by Data object arg over 
    its domain.

    @param arg:   Data object representing the function which is integrated.
    @param where: FunctionSpace in which the integral is calculated. 
                  If not present or C{None} an appropriate default is used.
    """
    if _testForZero(arg):
      return 0
    elif isinstance(arg,symbols.Symbol):
       return symbols.Integral_Symbol(arg,where)
    else:    
       if not where==None: arg=escript.Data(arg,where)
       if arg.getRank()==0:
         return arg.integrate()[0]
       else:
         return arg.integrate()

#=============================
#
# wrapper for various functions: if the argument has attribute the function name
# as an argument it calls the corresponding methods. Otherwise the corresponding
# numarray function is called.

# functions involving the underlying Domain:


# functions returning Data objects:

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
    if isinstance(arg,symbols.Symbol):
       return symbols.Transpose_Symbol(arg,axis=r)
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
    if isinstance(arg,symbols.Symbol):
       s=list(arg.getShape())        
       s=tuple(s[0:axis0]+s[axis0+1:axis1]+s[axis1+1:])
       return symbols.Trace_Symbol(arg,axis0=axis0,axis1=axis1)
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

def length(arg):
    """
    @param arg:
    """
    if isinstance(arg,escript.Data):
       if arg.isEmpty(): return escript.Data()
       if arg.getRank()==0:
          return abs(arg)
       elif arg.getRank()==1:
          out=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             out+=arg[i]**2
          return sqrt(out)
       elif arg.getRank()==2:
          out=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             for j in range(arg.getShape()[1]):
                out+=arg[i,j]**2
          return sqrt(out)
       elif arg.getRank()==3:
          out=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             for j in range(arg.getShape()[1]):
                for k in range(arg.getShape()[2]):
                   out+=arg[i,j,k]**2
          return sqrt(out)
       elif arg.getRank()==4:
          out=escript.Scalar(0,arg.getFunctionSpace())
          for i in range(arg.getShape()[0]):
             for j in range(arg.getShape()[1]):
                for k in range(arg.getShape()[2]):
                   for l in range(arg.getShape()[3]):
                      out+=arg[i,j,k,l]**2
          return sqrt(out)
       else:
          raise SystemError,"length is not been fully implemented yet"
          # return arg.length()
    elif isinstance(arg,float):
       return abs(arg)
    else:
       return sqrt((arg**2).sum())

def deviator(arg):
    """
    @param arg:
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
    @param arg0:
    @param arg1:
    """
    if isinstance(arg0,escript.Data):
       arg=arg0
    else:
       arg=arg1

    out=escript.Scalar(0,arg.getFunctionSpace())
    if arg.getRank()==0:
          return arg0*arg1
    elif arg.getRank()==1:
         out=escript.Scalar(0,arg.getFunctionSpace())
         for i in range(arg.getShape()[0]):
            out+=arg0[i]*arg1[i]
    elif arg.getRank()==2:
        out=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
           for j in range(arg.getShape()[1]):
              out+=arg0[i,j]*arg1[i,j]
    elif arg.getRank()==3:
        out=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
            for j in range(arg.getShape()[1]):
               for k in range(arg.getShape()[2]):
                  out+=arg0[i,j,k]*arg1[i,j,k]
    elif arg.getRank()==4:
        out=escript.Scalar(0,arg.getFunctionSpace())
        for i in range(arg.getShape()[0]):
           for j in range(arg.getShape()[1]):
              for k in range(arg.getShape()[2]):
                 for l in range(arg.getShape()[3]):
                    out+=arg0[i,j,k,l]*arg1[i,j,k,l]
    else:
          raise SystemError,"inner is not been implemented yet"
    return out

def tensormult(arg0,arg1):
    # check LinearPDE!!!!
    raise SystemError,"tensormult is not implemented yet!"

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
               # uses Data object slicing, plus Data * and += operators
               out[i]+=arg0[i,j]*arg1[j]
          return out
      elif arg0.getRank()==1 and arg1.getRank()==1:
          return inner(arg0,arg1)
      else:
          raise SystemError,"matrixmult is not fully implemented yet!"

#=========================================================
# reduction operations:
#=========================================================
def sum(arg):
    """
    @param arg:
    """
    return arg.sum()

def sup(arg):
    """
    @param arg:
    """
    if isinstance(arg,escript.Data):
       return arg.sup()
    elif isinstance(arg,float) or isinstance(arg,int):
       return arg
    else:
       return arg.max()

def inf(arg):
    """
    @param arg:
    """
    if isinstance(arg,escript.Data):
       return arg.inf()
    elif isinstance(arg,float) or isinstance(arg,int):
       return arg
    else:
       return arg.min()

def L2(arg):
    """
    Returns the L2-norm of the argument

    @param arg:
    """
    if isinstance(arg,escript.Data):
       return arg.L2()
    elif isinstance(arg,float) or isinstance(arg,int):
       return abs(arg)
    else:
       return numarry.sqrt(dot(arg,arg))

def Lsup(arg):
    """
    @param arg:
    """
    if isinstance(arg,escript.Data):
       return arg.Lsup()
    elif isinstance(arg,float) or isinstance(arg,int):
       return abs(arg)
    else:
       return numarray.abs(arg).max()

def dot(arg0,arg1):
    """
    @param arg0:
    @param arg1:
    """
    if isinstance(arg0,escript.Data):
       return arg0.dot(arg1)
    elif isinstance(arg1,escript.Data):
       return arg1.dot(arg0)
    else:
       return numarray.dot(arg0,arg1)

def kronecker(d):
   if hasattr(d,"getDim"):
      return numarray.identity(d.getDim())*1.
   else:
      return numarray.identity(d)*1.

def unit(i,d):
   """
   Return a unit vector of dimension d with nonzero index i.

   @param d: dimension
   @param i: index
   """
   e = numarray.zeros((d,),numarray.Float)
   e[i] = 1.0
   return e
