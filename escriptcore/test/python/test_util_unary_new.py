
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
test for util operations for unary operations without tagged data

:remark: use see `test_util`
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Joel Fenwick, joelfenwick@uq.edu.au"

import esys.escriptcore.utestselect as unittest
import numpy
import math
import cmath
from esys.escript import *
from test_util_base import Test_util_base, Test_util_values

haveLapack = hasFeature('lapack')

def getTagStr(f, p=None):
   if p is None:
        return "self.makeTagged(self.functionspace, "+f+"(r), "+f+"(r2))"
   else:
        f1=f.replace(p,'r')
        f2=f.replace(p,'r2')
        return "self.makeTagged(self.functionspace, ("+f1+"), ("+f2+"))"

class Test_util_unary_new(Test_util_values):
   """
   test for unary operations. 
   """
   def iterateops(self, ops, vals):
       """
       """
       for p in ops:
          o,c,z=p
          for v in vals:
            res=o(v)
            if isinstance(v,complex):
               ref=z(v)
            else:
               ref=c(v)
            self.assertTrue(isinstance(res,type(ref)),"wrong type of result for "+str(o))
            self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result for "+str(o))
            d=Data(v)
            res=o(d)
            self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result for data on "+str(o))

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_compare_complex_vs_real_data(self):
      # Compare results of unary ops provided by util and their python equivalents
      # Note that in some cases util calls these python versions so it is not a
      # guarantee that they are correct, we also compare with Data implementation
      # so we reduce the risk unless there is a fault in the system's underlying c library
      # Also note that these calls are only testing scalars
      ops=[(sin,math.sin,cmath.sin), (cos,math.cos,cmath.cos), (tan,math.tan,cmath.tan), (log,math.log,cmath.log), (log10, math.log10, cmath.log10), (Abs, abs, abs),
(acos,math.acos,cmath.acos), (acosh,math.acosh,cmath.acosh), (asin,math.asin,cmath.asin), (asinh, math.asinh,cmath.asinh),
(cosh, math.cosh, cmath.cosh), (exp, math.exp, cmath.exp), (sinh, math.sinh, cmath.sinh), (sqrt, math.sqrt, cmath.sqrt)]
      vals=[1+0j,-1+0j,1j, -1j, math.pi*1j,3+4j]
      self.iterateops(ops,vals)
      ops=[(atan,math.atan,cmath.atan)]
      vals=[1+0j,-1+0j, math.pi*1j,3+4j]
      self.iterateops(ops,vals)
      ops=[(atanh,math.atanh,cmath.atanh)]
      vals=[1j, -1j, math.pi*1j,3+4j]
      self.iterateops(ops,vals)
      # test for zero values for those functions which can take it
      vals=[0j]
      ops=[(sin,math.sin,cmath.sin), (cos,math.cos,cmath.cos), (tan,math.tan,cmath.tan),
           (Abs, abs, abs),
           (acos,math.acos,cmath.acos), (acosh,math.acosh,cmath.acosh),
           (asin,math.asin,cmath.asin), (asinh, math.asinh,cmath.asinh),
           (cosh, math.cosh, cmath.cosh), (exp, math.exp, cmath.exp), 
           (sinh, math.sinh, cmath.sinh), (sqrt, math.sqrt, cmath.sqrt),
           (atan,math.atan,cmath.atan), (atanh,math.atanh,cmath.atanh)]
      self.iterateops(ops,vals)      

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_abs_combined(self):
        supportcplx=True
        opstring="abs(a)"
        misccheck="isinstance(res,float) or (isinstance(res, numpy.ndarray) and res.dtype.kind=='f') or (isinstance(res, Data) and not res.isComplex())"
        oraclecheck="numpy.abs(ref)"
        opname="abs"
        update1=getTagStr("numpy.abs")    
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1)  
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_acos_combined(self):
        supportcplx=True
        opstring="acos(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.arccos(ref)"
        opname="acos"
        update1=getTagStr("numpy.arccos")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_acosh_combined(self):
        supportcplx=True
        opstring="acosh(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.arccosh(ref)"
        opname="acosh"
        update1=getTagStr("numpy.arccosh")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, input_trans=lambda x: x+2)  
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_antihermitian_combined(self):
        supportcplx=True
        opstring="antihermitian(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="(ref-numpy.conjugate(transpose(ref)))/2"
        opname="antihermitian"
        update1=getTagStr("(r2-numpy.conjugate(transpose(r2)))/2",p='r2')
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, expect_raise_on_ranks=(0,1,3), expected_exceptions=(ValueError,))         
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_antisymmetric_combined(self):
        supportcplx=True
        opstring="antisymmetric(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="(ref-transpose(ref))/2"
        opname="antisymmetric"
        update1=getTagStr("(r2-transpose(r2))/2",p='r2')
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, expect_raise_on_ranks=(0,1,3), expected_exceptions=(ValueError,))            
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_asin_combined(self):
        supportcplx=True
        opstring="asin(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.arcsin(ref)"
        opname="asin"
        update1=getTagStr("numpy.arcsin")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_asinh_combined(self):
        supportcplx=True
        opstring="asinh(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.arcsinh(ref)"
        opname="asinh"
        update1=getTagStr("numpy.arcsinh")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)   
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_atan_combined(self):
        supportcplx=True
        opstring="atan(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.arctan(ref)"
        opname="atan"
        update1=getTagStr("numpy.arctan")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)   
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_atanh_combined(self):
        supportcplx=True
        opstring="atanh(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.arctanh(ref)"
        opname="atanh"
        update1=getTagStr("numpy.arctanh")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1) 
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_clip_combined(self):
        supportcplx=False
        opstring="clip(a,minval=-0.5, maxval=0.5)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.clip(ref, -0.5, 0.5)"
        opname="clip"
        update1=getTagStr("numpy.clip(r2, -0.5, 0.5)", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)          
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_cos_combined(self):
        supportcplx=True
        opstring="cos(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.cos(ref)"
        opname="cos"
        update1=getTagStr("numpy.cos")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)  
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_cosh_combined(self):
        supportcplx=True
        opstring="cosh(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.cosh(ref)"
        opname="cosh"
        update1=getTagStr("numpy.cosh")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)  
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_eigenvalues_combined(self):
        supportcplx=True
        opstring="eigenvalues(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.linalg.eigvals(ref) if numpy.ndim(ref)==2 else 0."
        opname="eigenvalues"
        update1=getTagStr("numpy.sort(numpy.linalg.eigvals(r2)) if numpy.ndim(r2)==2 else 0.", p='r2')
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1,
                                                 expect_raise_on_ranks=(0,1,3,4), no_scalars=True, 
                                                 input_trans=lambda x: symmetric(x) if getRank(x)==2 else x)
        # We aren't testing the scalar case
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_eigenvalues_and_eigenvectors_combined(self):
        supportcplx=False
        opstring="eigenvalues_and_eigenvectors(a)"
        misccheck="(Lsup(matrix_mult(a,res[1][:,0])-res[0][0]*res[1][:,0])<=0.001*Lsup(res[0])) and Lsup(matrix_mult(a,res[1][:,1])-res[0][1]*res[1][:,1])<=0.001*Lsup(res[0])"
        oraclecheck=None
        opname="eigenvalues_and_eigenvectors"
        update1=None    # note missing update1 step, hoping this won't be a problem
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1,
                                                 expect_raise_on_ranks=(0,1,3,4), no_scalars=True,
                                                 input_trans=lambda x: symmetric(x) if getRank(x)==2 else x)
        # We aren't testing the scalar case        
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_exp_combined(self):
        supportcplx=True
        opstring="exp(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.exp(ref)"
        opname="exp"
        update1=getTagStr("numpy.exp")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_hermetitian_combined(self):
        supportcplx=True
        opstring="hermitian(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="(ref+numpy.conjugate(transpose(ref)))/2"
        opname="hermitian"
        update1=getTagStr("(r2+numpy.conjugate(transpose(r2)))/2", p='r2')
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
                                            expect_raise_on_ranks=(0,1,3), expected_exceptions=(ValueError,))
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inverse_combined(self):
        supportcplx=False
        opstring="inverse(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.linalg.inv(ref)"
        opname="inverse"
        update1=getTagStr("numpy.linalg.inv")
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
                                            minrank=2, maxrank=2, no_scalars=True)
        # escript's inverse also supports scalars so need to check them separately
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_phase_combined(self):
        supportcplx=True
        opstring="phase(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.angle(ref)"
        opname="phase"
        update1=getTagStr("numpy.angle")
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
                                            minrank=2, maxrank=2, no_scalars=True)        
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sign_combined(self):
        supportcplx=False
        opstring="sign(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.sign(ref)"
        opname="sign"
        update1=getTagStr("numpy.sign")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)          
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sin_combined(self):
        supportcplx=True
        opstring="sin(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.sin(ref)"
        opname="sin"
        update1=getTagStr("numpy.sin")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)  
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sinh_combined(self):
        supportcplx=True
        opstring="sinh(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.sinh(ref)"
        opname="sinh"
        update1=getTagStr("numpy.sinh")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)  
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sqrt_combined(self):
        supportcplx=True
        opstring="sqrt(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.sqrt(ref)"
        opname="sqrt"
        update1=getTagStr("numpy.sqrt")  
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, input_trans=lambda x: numpy.abs(x) if type(x) is numpy.ndarray and x.dtype.kind=='f' else abs(x) if type(x) is Data and not x.isComplex() else x)          
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_symmetric_combined(self):
        supportcplx=True
        opstring="symmetric(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="(ref+transpose(ref))/2"
        opname="symmetric"
        update1=getTagStr("(r2+transpose(r2))/2", p='r2')
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
                                            expect_raise_on_ranks=(0,1,3), expected_exceptions=(ValueError,))         
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_swapaxes_combined(self):
        supportcplx=True
        opstring="swap_axes(a,axis0=0, axis1=getRank(a)-1)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.swapaxes(ref, 0,ref.ndim-1)"
        opname="swap_axes"
        update1=getTagStr("numpy.swapaxes(r2, 0, r2.ndim-1)", p='r2')
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
                                            minrank=2, maxrank=4,no_scalars=True)                
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_tan_combined(self):
        supportcplx=True
        opstring="tan(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.tan(ref)"
        opname="tan"
        update1=getTagStr("numpy.tan")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1) 
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_tanh_combined(self):
        supportcplx=True
        opstring="tanh(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.tanh(ref)"
        opname="tanh"
        update1=getTagStr("numpy.tanh") 
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_trace_combined(self):
        supportcplx=True
        opstring="trace(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.trace(ref)"
        opname="trace"
        update1=getTagStr("numpy.trace(r2) if numpy.ndim(r2)>=2 else None", p='r2')
        update2=None
        # We could also check to see if it throws when given a scalar but we don't
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
                                      no_scalars=True, minrank=2,
                                           expect_raise_on_ranks=(0,1),expected_exceptions=(ValueError,))
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_transpose_combined(self):
        supportcplx=True
        opstring="transpose(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.transpose(ref, axes=list(range(len(ref.shape)//2,len(ref.shape)))+list(range(0,len(ref.shape)//2)))"
        opname="tranpose"
        update1=getTagStr("numpy.transpose(r2, axes=list(range(len(r2.shape)//2,len(r2.shape)))+list(range(0,len(r2.shape)//2)))", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)   
        opname="transpose-offset-1"
        opstring="transpose(a,axis_offset=1)"
        oraclecheck="numpy.transpose(ref, axes=list(range(1,len(ref.shape)))+list(range(0,1)))"
        update1=getTagStr("numpy.transpose(r2, axes=list(range(1,len(r2.shape)))+list(range(0,1)))", p='r2')
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1,
                                      no_scalars=True)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_length_combined(self):
        supportcplx=True
        opstring="length(a)"
        misccheck="isinstance(res,float) or (isinstance(res,Data) and not res.isComplex())"
        oraclecheck="numpy.linalg.norm(ref)"
        opname="length"
        update1=getTagStr("numpy.linalg.norm")
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)         
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_log_combined(self):
        supportcplx=True
        opstring="log(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.log(ref)"
        opname="log"
        update1=getTagStr("numpy.log")
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, 
            input_trans=lambda x: numpy.abs(x)+1 if type(x) is numpy.ndarray and x.dtype.kind=='f' else abs(x)+1 if type(x) is Data and not x.isComplex() else x+1)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_log10_combined(self):
        supportcplx=True
        opstring="log10(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.log10(ref)"
        opname="log10"
        update1=getTagStr("numpy.log10") 
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1,  
            input_trans=lambda x: numpy.abs(x)+1 if type(x) is numpy.ndarray and x.dtype.kind=='f' else abs(x)+1 if type(x) is Data and not x.isComplex() else x+1)  
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_maxval_combined(self):
        supportcplx=False
        opstring="maxval(a)"
        misccheck="isinstance(res,float) or (isinstance(a,numpy.ndarray) and res.dtype.kind=='f') or (isinstance(a, Data) and not res.isComplex())"
        oraclecheck="numpy.max(ref)"
        opname="maxval"
        update1=getTagStr("numpy.max") 
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1)  
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_minval_combined(self):
        supportcplx=False
        opstring="minval(a)"
        misccheck="isinstance(res,float) or (isinstance(a,numpy.ndarray) and res.dtype.kind=='f') or (isinstance(a, Data) and not res.isComplex())"
        oraclecheck="numpy.min(ref)"
        opname="minval"
        update1=getTagStr("numpy.min")    
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1) 
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_whereNegative_combined(self):
        supportcplx=False
        opstring="whereNegative(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.where(ref<0, numpy.ones(ref.shape), numpy.zeros(ref.shape))"
        opname="whereNegative"
        update1=getTagStr("numpy.where(r2<0, numpy.ones(r2.shape), numpy.zeros(r2.shape))", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)     
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_whereNonNegative_combined(self):
        supportcplx=False
        opstring="whereNonNegative(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.where(ref>=0, numpy.ones(ref.shape), numpy.zeros(ref.shape))"
        opname="whereNonNegative"
        update1=getTagStr("numpy.where(r2>=0, numpy.ones(r2.shape), numpy.zeros(r2.shape))", p='r2')    
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)    
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_whereNonPositive_combined(self):
        supportcplx=False
        opstring="whereNonPositive(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.where(ref<=0, numpy.ones(ref.shape), numpy.zeros(ref.shape))"
        opname="whereNonPositive"
        update1=getTagStr("numpy.where(r2<=0, numpy.ones(r2.shape), numpy.zeros(r2.shape))", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_whereNonZero_combined(self):
        supportcplx=True
        opstring="whereNonZero(a)"
        misccheck="isinstance(res,float) or (isinstance(res, numpy.ndarray) and res.dtype.kind=='f') or (isinstance(res, Data) and not res.isComplex())"
        oraclecheck="numpy.where(ref!=0, numpy.ones(ref.shape), numpy.zeros(ref.shape))"
        opname="whereNonZero"
        update1=getTagStr("numpy.where(r2!=0, numpy.ones(r2.shape), numpy.zeros(r2.shape))", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)        
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_wherePositive_combined(self):
        supportcplx=False
        opstring="wherePositive(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.where(ref>0, numpy.ones(ref.shape), numpy.zeros(ref.shape))"
        opname="wherePositive"
        update1=getTagStr("numpy.where(r2>0, numpy.ones(r2.shape), numpy.zeros(r2.shape))", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1) 
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_whereZero_combined(self):
        supportcplx=True
        opstring="whereZero(a)"
        misccheck="isinstance(res,float) or (isinstance(res, numpy.ndarray) and res.dtype.kind=='f') or (isinstance(res, Data) and not res.isComplex())"
        oraclecheck="numpy.where(ref==0, numpy.ones(ref.shape), numpy.zeros(ref.shape))"
        opname="whereZero"
        update1=getTagStr("numpy.where(r2==0, numpy.ones(r2.shape), numpy.zeros(r2.shape))", p='r2')
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
