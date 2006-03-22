// $Id$
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#if !defined escript_UnaryOp_20040315_H
#define escript_UnaryOp_20040315_H

#include "DataArrayView.h"
#include "DataConstant.h"
#include "DataTagged.h"
#include "DataExpanded.h"

namespace escript {

/**
   \brief
   Perform the given unary operation on each data point of the given Data object.
   Called by Data::unaryOp.
   Calls DataArrayView::unaryOp.
   For DataExpanded objects, operation is done in parallel.
   \param data Input/Output - The data.
   \param operation Input - The operation to perform.
*/

template <class UnaryFunction>
inline
void
unaryOp(DataExpanded& data,
        UnaryFunction operation)
{
  int i,j;
  DataArrayView::ValueType::size_type numDPPSample=data.getNumDPPSample();
  DataArrayView::ValueType::size_type numSamples=data.getNumSamples();
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numSamples;i++) {
    for (j=0;j<numDPPSample;j++) {
      data.getPointDataView().unaryOp(data.getPointOffset(i,j),operation);
    }
  }
}

template <class UnaryFunction>
inline
void
unaryOp(DataTagged& data,
        UnaryFunction operation)
{
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataArrayView& dataView=data.getPointDataView();
  for (i=lookup.begin();i!=lookupEnd;i++) {
    dataView.unaryOp(i->second,operation);
  }
  // perform the operation on the default value
  data.getDefaultValue().unaryOp(operation);
}

template <class UnaryFunction>
inline
void
unaryOp(DataConstant& data,
        UnaryFunction operation)
{
  data.getPointDataView().unaryOp(operation);
}

} // end of namespace
#endif
