// $Id$
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
                                                                           
#if !defined escript_DataAlgorithm_20040714_H
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
   Return the length.

   Description:
   Return the length.
*/
struct Length : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return std::sqrt(std::pow(x,2)+std::pow(y,2));
  }
};

/**
   \brief
   Return the trace.

   Description:
   Return the trace.
*/
struct Trace : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return x+y;
  }
};

/**
   \brief
   Adapt algorithms so they may be used by Data.

   Description:
   Adapt algorithms so they may be used by Data. The functor 
   maintains state, the currentValue returned by the operation,
   and the initial value.
*/
template <class BinaryFunction>
class DataAlgorithmAdapter {
 public:
    DataAlgorithmAdapter(double initialValue):
      m_initialValue(initialValue),
      m_currentValue(initialValue)
    {}
    inline void operator()(double value)
    {
      m_currentValue=operation(m_currentValue,value);
      return;
    }
    inline void resetResult()
    {
      m_currentValue=m_initialValue;
    }
    inline double getResult() const
    {
      return m_currentValue;
    }
 private:
    //
    // the initial operation value
    double m_initialValue;
    //
    // the current operation value
    double m_currentValue;
    //
    // The operation to perform
    BinaryFunction operation;
};

/**
   \brief
   Perform the given operation upon all Data elements and return a single
   result.

   Description:
   Perform the given operation upon all Data elements and return a single
   result.
*/
template <class UnaryFunction>
inline
double
algorithm(DataExpanded& data,
          UnaryFunction operation)
{
  int i,j;
  DataArrayView::ValueType::size_type numDPPSample=data.getNumDPPSample();
  DataArrayView::ValueType::size_type numSamples=data.getNumSamples();
  double resultLocal=0;
#pragma omp parallel private(resultLocal)
  {
#pragma omp for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	resultLocal=data.getPointDataView().reductionOp(data.getPointOffset(i,j), operation);
#pragma omp critical (algorithm)
	operation(resultLocal);
      }
    }
  }
  return operation.getResult();
}

template <class UnaryFunction>
inline
double
algorithm(DataTagged& data,
          UnaryFunction operation)
{
  //
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataArrayView& dataView=data.getPointDataView();
  for (i=lookup.begin();i!=lookupEnd;i++) {
    operation(dataView.reductionOp(i->second,operation));
  }
  //
  // finally perform the operation on the default value
  operation(data.getDefaultValue().reductionOp(operation));
  return operation.getResult();
}

template <class UnaryFunction>
inline
double
algorithm(DataConstant& data,
          UnaryFunction operation)
{
  return data.getPointDataView().reductionOp(operation);
}

/**
   \brief
   Perform the given data point reduction operation on all data points
   in data, storing results in corresponding elements of result.

   Objects data and result must be of the same type, and have the same number
   of data points, but where data has data points of rank n, result must have
   data points of rank 0.

   Calls DataArrayView::dp_algorithm.
*/
template <class UnaryFunction>
inline
void
dp_algorithm(DataExpanded& data,
             DataExpanded& result,
             UnaryFunction operation)
{
  //
  // perform the operation on each data value
  // and assign this to the corresponding element in result
  int i,j;
  DataArrayView::ValueType::size_type numDPPSample=data.getNumDPPSample();
  DataArrayView::ValueType::size_type numSamples=data.getNumSamples();
  {
#pragma omp for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
#pragma omp critical (dp_algorithm)
        result.getPointDataView().getData(data.getPointOffset(i,j)) =
          data.getPointDataView().dp_reductionOp(data.getPointOffset(i,j),operation);
      }
    }
  }
}

template <class UnaryFunction>
inline
void
dp_algorithm(DataTagged& data,
             DataTagged& result,
             UnaryFunction operation)
{
  //
  // perform the operation on each tagged data value
  // and assign this to the corresponding element in result
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  for (i=lookup.begin();i!=lookupEnd;i++) {
    result.getPointDataView().getData(i->second) =
      data.getPointDataView().dp_reductionOp(i->second,operation);
  }
  //
  // finally perform the operation on the default data value
  // and assign this to the default element in result
  result.getPointDataView().getData(0) =
    data.getDefaultValue().dp_reductionOp(operation);
}

template <class UnaryFunction>
inline
void
dp_algorithm(DataConstant& data,
             DataConstant& result,
             UnaryFunction operation)
{
  //
  // perform the operation on the default data value
  // and assign this to the default element in result
  result.getPointDataView().getData(0) =
    data.getPointDataView().dp_reductionOp(operation);
}

} // end of namespace
#endif
