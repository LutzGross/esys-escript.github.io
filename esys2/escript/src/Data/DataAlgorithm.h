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
                                                                           
#if !defined  escript_DataAlgorithm_20040714_H
#define escript_DataAlgorithm_20040714_H

#include "escript/Data/DataExpanded.h"
#include "escript/Data/DataTagged.h"
#include "escript/Data/DataConstant.h"
#include "escript/Data/DataArrayView.h"

#include <iostream>
#include <algorithm>
#include <math.h>
#include <limits>

namespace escript {
/**
   \brief
   Return the maximum value.

   Description:
   Return the maximum value.
*/
struct FMax : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return std::max(x,y);
  }
};

/**
   \brief
   Return the minimum value.

   Description:
   Return the minimum value.
*/
struct FMin : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return std::min(x,y);
  }
};

/**
   \brief
   Return the absolute maximum value.

   Description:
   Return the absolute maximum value.
*/
struct AbsMax : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return std::max(fabs(x),fabs(y));
  }
};
/**
   \brief
   Adapt algorithms so they may be used by Data.

   Description:
   Adapt algorithms so they may be used by Data. The functor 
   maintains state, ie the currentValue retuned by the operation.
*/
template <class BinaryFunction>
class DataAlgorithmAdapter {
 public:
    DataAlgorithmAdapter(double initialValue):
      m_currentValue(initialValue)
    {}
    inline void operator()(double value)
    {
      m_currentValue=operation(m_currentValue,value);
      return;
    }
    inline double getResult() const
    {
      return m_currentValue;
    }
 private:
    //
    // the current maximum value
    double m_currentValue;
    //
    // The operation to perform
    BinaryFunction operation;
};

/**
   \brief
   Perform the given operation upon all DataElements and return a single
   result.

   Description:
   Perform the given operation upon all DataElements and return a single
   result.
*/
template <class UnaryFunction>
inline double algorithm(DataExpanded& data, UnaryFunction operation)
{
  int i,j;
  DataArrayView::ValueType::size_type numDPPSample=data.getNumDPPSample();
  DataArrayView::ValueType::size_type numSamples=data.getNumSamples();
  double resultLocal;
  #pragma omp parallel private(resultLocal)
  {
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;++i) {
      for (j=0;j<numDPPSample;++j) {
	resultLocal=data.getPointDataView().algorithm(data.getPointOffset(i,j),
						      operation);
	#pragma omp critical (algorithm)
	operation(resultLocal);
      }
    }
  }
  return operation.getResult();
}

template <class UnaryFunction>
inline double algorithm(DataTagged& data, UnaryFunction operation)
{
  //
  // perform the operation on each tagged value including the default
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataArrayView& dataView=data.getPointDataView();
  for (i=lookup.begin();i!=lookupEnd;++i) {
    operation(dataView.algorithm(i->second,operation));
  }
  //
  // finally perform the operation on the default value
  operation(data.getDefaultValue().algorithm(operation));
  return operation.getResult();
}

template <class UnaryFunction>
inline double algorithm(DataConstant& data, UnaryFunction operation)
{
  return data.getPointDataView().algorithm(operation);
}


} // end of namespace
#endif
