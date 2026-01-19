
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
test for non-overloaded binary operations

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
from esys.escript import *
from test_util_base import Test_util_values




        
    

class Test_util_binary_new(Test_util_values):
    
   def generate_indices(self, shape):
        res=[0]*len(shape)
        l=len(shape)
        done=False
        while not done:
            yield tuple(res)
            res[0]+=1
            for i in range(l-1):
                if res[i]>=shape[i]:
                    res[i]=0
                    res[i+1]+=1
                else:
                    break
            # now we check the last digit
            if res[l-1]>=shape[l-1]:
                done=True
        

   def subst_outer(self, a, b):
        if isinstance(a,float) or isinstance(a, complex):
            a=(a,)
        if isinstance(b,float) or isinstance(b, complex):
            b=(b,)
        sa=getShape(a)
        sb=getShape(b)
        a=numpy.array(a)
        b=numpy.array(b)
        targettype=a.dtype if a.dtype.kind=='c' else b.dtype
        if sa==():
            if sb==():
                return a*b
            resshape=sb
            res=numpy.zeros(resshape, dtype=targettype)
            for xb in self.generate_indices(sb):
                res[xb] = a*b[xb]
            return res
        elif sb==():
            resshape=sa
            res=numpy.zeros(resshape, dtype=targettype)
            for xa in self.generate_indices(sa):
                res[xa] = a[xa]*b
            return res
        else:
            resshape=sa+sb
            res=numpy.zeros(resshape, dtype=targettype)
        for xa in self.generate_indices(sa):
            for xb in self.generate_indices(sb):
                res[xa+xb] = a[xa]*b[xb]
        return res    
    
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_combined(self):
       opstring='a+b'
       misccheck='isinstance(res, Data) if isinstance(a, Data) or isinstance(b, Data) else True' # doesn't cover all cases;
       oraclecheck="refa+refb"
       opname="+ operator"
       noshapemismatch=True
       permitscalarmismatch=True
       permit_array_op_data=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, permit_array_op_data=permit_array_op_data)   
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_combined(self):
       opstring='a/b'
       misccheck='isinstance(res, Data) if isinstance(a, Data) or isinstance(b, Data) else True' # doesn't cover all cases;
       oraclecheck="refa/refb"
       opname="/ operator"
       no_second_arg_zeros=True
       noshapemismatch=True
       permitscalarmismatch=True
       permit_array_op_data=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, permit_array_op_data=permit_array_op_data, no_second_arg_zeros=no_second_arg_zeros)     
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_mul_combined(self):
       opstring='a*b'
       misccheck='isinstance(res, Data) if isinstance(a, Data) or isinstance(b, Data) else True' # doesn't cover all cases;
       oraclecheck="refa*refb"
       opname="* operator"
       noshapemismatch=True
       permitscalarmismatch=True
       permit_array_op_data=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, permit_array_op_data=permit_array_op_data)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_pow_combined(self):
       opstring='a**b'
       misccheck='isinstance(res, Data) if isinstance(a, Data) or isinstance(b, Data) else True' # doesn't cover all cases;
       oraclecheck="refa**refb"
       opname="** operator"
       noshapemismatch=True
       permitscalarmismatch=True
       permit_array_op_data=False
       no_first_arg_negative=True
       no_first_arg_zeros=False
       second_large_args=False
       first_large_args=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, permit_array_op_data=permit_array_op_data, no_first_arg_negative=no_first_arg_negative, no_first_arg_zeros=no_first_arg_zeros, 
                                                 second_large_args=second_large_args, first_large_args=first_large_args)          
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sub_combined(self):
       opstring='a-b'
       misccheck='isinstance(res, Data) if isinstance(a, Data) or isinstance(b, Data) else True' # doesn't cover all cases;
       oraclecheck="refa-refb"
       opname="- operator"
       noshapemismatch=True
       permitscalarmismatch=True
       permit_array_op_data=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, permit_array_op_data=permit_array_op_data)          
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inner_combined(self):
       opstring='inner(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.tensordot(refa, refb, axes=getRank(refa))"
       opname="inner"
       noshapemismatch=True
       permitscalarmismatch=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_outer_combined(self):
       opstring='outer(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="self.subst_outer(refa,refb)"
       opname="outer"
       noshapemismatch=True
       permitscalarmismatch=True
       capcombinedrank=True
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, cap_combined_rank=capcombinedrank)           
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_matrix_minimum_combined(self):
       opstring='minimum(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.minimum(refa,refb)"
       opname="minimum"
       noshapemismatch=True
       permitscalarmismatch=True
       cplx=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, support_cplx=cplx)         
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_matrix_maximum_combined(self):
       opstring='maximum(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.maximum(refa,refb)"
       opname="maximum"
       noshapemismatch=True
       permitscalarmismatch=True
       cplx=False
       self.generate_binary_operation_test_batch(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch, support_cplx=cplx)         
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_matrix_mult_combined(self):
       opstring='matrix_mult(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.dot(refa,refb)"
       opname="matrix_mult"
       aranks=(2,)
       self.generate_binary_matrixlike_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, aranks=aranks)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_transpose_matrix_mult_combined(self):
       opstring='transposed_matrix_mult(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.dot(numpy.transpose(refa),refb)"
       opname="transposed_matrix_mult"
       aranks=(2,)
       self.generate_binary_matrixlike_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, aranks=aranks)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_matrix_transposed_mult_combined(self):
       opstring='matrix_transposed_mult(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.dot(refa,numpy.transpose(refb))"
       opname="matrix_transposed_mult"
       aranks=(2,)
       self.generate_binary_matrixlike_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, aranks=aranks)       
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_tensor_mult_combined(self):
       opstring='tensor_mult(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.dot(refa,refb) if getRank(refa)==2 else numpy.tensordot(refa,refb)"
       opname="tensor_mult"
       aranks=(2,4)
       self.generate_binary_matrixlike_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, aranks=aranks)
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_transposed_tensor_mult_combined(self):
       opstring='transposed_tensor_mult(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.dot(transpose(refa),refb) if getRank(refa)==2 else numpy.tensordot(transpose(refa),refb)"
       opname="transposed_tensor_mult"
       aranks=(2,4)
       self.generate_binary_matrixlike_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, aranks=aranks) 
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_tensor_transposed_mult_combined(self):
       opstring='tensor_transposed_mult(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.dot(refa,transpose(refb)) if getRank(refa)==2 else numpy.tensordot(refa,transpose(refb))"
       opname="tensor_tranposed_mult"
       aranks=(2,4)
       self.generate_binary_matrixlike_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, aranks=aranks)         
