/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "RTOpPack_ROpGetSubVector.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetSubVector, apply_op_1, Scalar )
{

  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(2*n, ST::zero());
  for (index_type k = 0; k < n; ++k)
    sv(k) = ST::random();

  RTOpPack::ROpGetSubVector<Scalar> op(0, 2*n-1);

  RCP<RTOpPack::ReductTarget> reductObj = rtopt(op).reduct_obj_create();

  rtopt(op).apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    reductObj.ptr()
    );

  const ConstSubVectorView<Scalar> reduct_sv = op(*reductObj);

  if (verbose) {
    dumpSubVectorView(sv, "sv", out);
    dumpSubVectorView(reduct_sv, "reduct_sv", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(sv),
    constSubVectorViewAsArray(reduct_sv) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetSubVector, apply_op_1 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetSubVector, apply_op_2, Scalar )
{

  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv0 = newSubVectorView<Scalar>(n, ST::zero());
  for (index_type k = 0; k < n; ++k)
    sv0(k) = ST::random();

  SubVectorView<Scalar> sv1 = newSubVectorView<Scalar>(n, ST::zero());
  sv1.setGlobalOffset(n);
  for (index_type k = 0; k < n; ++k)
    sv1(k) = ST::random();

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(2*n, ST::zero());
  for (index_type k = 0; k < n; ++k)
    sv[k] = sv0(k);
  for (index_type k = 0; k < n; ++k)
    sv[k+n] = sv1(k);

  RTOpPack::ROpGetSubVector<Scalar> op(0, 2*n-1);

  RCP<RTOpPack::ReductTarget> reductObj = rtopt(op).reduct_obj_create();

  rtopt(op).apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    reductObj.ptr()
    );

  const ConstSubVectorView<Scalar> reduct_sv = op(*reductObj);

  if (verbose) {
    dumpSubVectorView(sv0, "sv0", out);
    dumpSubVectorView(sv1, "sv1", out);
    dumpSubVectorView(sv, "sv", out);
    dumpSubVectorView(reduct_sv, "reduct_sv", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(sv),
    constSubVectorViewAsArray(reduct_sv) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetSubVector, apply_op_2 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetSubVector, reductObjState, Scalar )
{

  using Teuchos::Array;
  using Teuchos::outArg;
  typedef ScalarTraits<Scalar> ST;
  typedef RTOpPack::PrimitiveTypeTraits<Scalar,Scalar> PTT;

  SubVectorView<Scalar> sv0 = newSubVectorView<Scalar>(n, ST::zero());
  for (index_type k = 0; k < n; ++k)
    sv0(k) = ST::random();

  RTOpPack::ROpGetSubVector<Scalar> op(0, n-1);
  
  int num_values = -1, num_indexes = -1, num_chars = -1;
  rtopt(op).get_reduct_type_num_entries(
    outArg(num_values), outArg(num_indexes), outArg(num_chars) );

  Array<typename PTT::primitiveType> value_data(num_values);
  Array<index_type> index_data(num_indexes);
  Array<char> char_data(num_chars);

  RCP<ReductTarget> reduct_obj_0 = rtopt(op).reduct_obj_create();

  rtopt(op).apply_op( tuple(ConstSubVectorView<Scalar>(sv0))(),
    null, outArg(*reduct_obj_0) );

  rtopt(op).extract_reduct_obj_state(
    *reduct_obj_0, value_data(), index_data(), char_data() );

  RCP<ReductTarget> reduct_obj_1 = rtopt(op).reduct_obj_create();

  rtopt(op).load_reduct_obj_state(
    value_data(), index_data(), char_data(),
    reduct_obj_1()
    );
  
  const ConstSubVectorView<Scalar> sv1 = op(*reduct_obj_1);

  if (verbose) {
    dumpSubVectorView(sv0, "sv0", out);
    dumpSubVectorView(sv1, "sv1", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(sv0),
    constSubVectorViewAsArray(sv1) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetSubVector, reductObjState )


} // namespace
