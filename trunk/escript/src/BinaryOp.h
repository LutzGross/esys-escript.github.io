/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined  escript_BinaryOp_20040315_H
#define escript_BinaryOp_20040315_H

#include "DataException.h"
#include "DataArrayView.h"
#include "DataConstant.h"
#include "DataExpanded.h"

#include <boost/scoped_ptr.hpp>

#include <iostream>
#include <functional>
#include <string>

namespace escript {
/**
   \brief
   Perform the given binary operation.
   \param left Input/Output - The left hand side.
   \param right Input - The right hand side.
   \param operation Input - The operation to perform.
*/
template <class BinaryFunction>
inline void binaryOp(DataTagged& left, const DataConstant& right, 
		     BinaryFunction operation)
{
  //
  // perform the operation on each tagged value including the default
  binaryOp(left,right.getPointDataView(),operation);
}

template <class BinaryFunction>
inline void binaryOp(DataTagged& left, const DataArrayView& right, 
		     BinaryFunction operation)
{
  //
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataArrayView& leftView=left.getPointDataView();
  if (right.getRank()==0) {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      leftView.binaryOp(i->second,right(),operation);
    }
  } else {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      leftView.binaryOp(i->second,right,right.getOffset(),operation);
    }
  }
  //
  // finally perform the operation on the default value
  left.getDefaultValue().binaryOp(right,operation);
}

template <class BinaryFunction>
inline void binaryOp(DataTagged& left, const DataTagged& right, 
		     BinaryFunction operation)
{
  //
  // Add the right hand tag keys which can't currently be found on the left
  const DataTagged::DataMapType& rightLookup=right.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator rightLookupEnd=rightLookup.end();
  for (i=rightLookup.begin();i!=rightLookupEnd;i++) {
    //
    // If the left does not already have a value assigned to this tag,
    // add the right hand tag to the left hand tag list and assign
    // the right's default value.
    if (!left.isCurrentTag(i->first)) {
      left.addTaggedValue(i->first,right.getDefaultValue());
    }
  }
  //
  // Perform the operation. Any tags originally in the left which don't exist for 
  // the right hand side will use the right's default value as the right operand
  const DataTagged::DataMapType& leftLookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator leftLookupEnd=leftLookup.end();
  for (i=leftLookup.begin();i!=leftLookupEnd;i++) {
    left.getDataPointByTag(i->first).binaryOp(right.getDataPointByTag(i->first),operation);
  }
  //
  // finally perform the operation on the default value
  left.getDefaultValue().binaryOp(right.getDefaultValue(),operation);
}

template <class BinaryFunction>
inline void binaryOp(DataConstant& left, const DataConstant& right, 
		     BinaryFunction operation)
{
  //
  // As DataConstant there is only one point data
  binaryOp(left,right.getPointDataView(),operation);
}

template <class BinaryFunction>
inline void binaryOp(DataConstant& left, const DataArrayView& right, 
		     BinaryFunction operation)
{
  //
  // perform an operand check, this will throw on error
  if (right.getRank()==0) {
    //
    // special case of applying a single value to the entire array
    left.getPointDataView().binaryOp(right(),operation);
  } else {
    left.getPointDataView().binaryOp(right,operation);
  }
}

template <class BinaryFunction>
inline void binaryOp(DataExpanded& left, const DataAbstract& right, 
		     BinaryFunction operation)
{
  int i,j;
  DataArrayView::ValueType::size_type numDPPSample=left.getNumDPPSample();
  DataArrayView::ValueType::size_type numSamples=left.getNumSamples();
  if (right.getPointDataView().getRank()==0) {
    //
    // This will call the double version of binaryOp
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	left.getPointDataView().binaryOp(left.getPointOffset(i,j),
					 right.getPointDataView().getData(right.getPointOffset(i,j)),
                                         operation);
      }
    }
  } else {
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	left.getPointDataView().binaryOp(left.getPointOffset(i,j),
					 right.getPointDataView(),
					 right.getPointOffset(i,j),
					 operation);
      }
    }
  }
}

template <class BinaryFunction>
inline void binaryOp(DataExpanded& left, const DataArrayView& right, 
		     BinaryFunction operation)
{
  int i,j;
  DataArrayView::ValueType::size_type numDPPSample=left.getNumDPPSample();
  DataArrayView::ValueType::size_type numSamples=left.getNumSamples();
  if (right.getRank()==0) {
    //
    // This will call the double version of binaryOp
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
       for (j=0;j<numDPPSample;j++) {
	left.getPointDataView().binaryOp(left.getPointOffset(i,j),
					 right(),
                                         operation);
      }
    }
  } else {
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	left.getPointDataView().binaryOp(left.getPointOffset(i,j),
					 right,
                                         0,
					 operation);
      }
    }
  }
}

} // end of namespace

#endif
