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
#include <vector>

namespace escript {

/**
   \brief
   Adapt binary algorithms so they may be used in DataArrayView reduction operations.

   Description:
   This functor adapts the given BinaryFunction operation by starting with the
   given inital value applying this operation to successive values, storing the
   rolling result in m_currentValue - which can be accessed or reset by getResult
   and resetResult respectively.
*/
template <class BinaryFunction>
class DataAlgorithmAdapter {
  public:
    DataAlgorithmAdapter(double initialValue):
      m_initialValue(initialValue),
      m_currentValue(initialValue)
    {
    }
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
   Return the maximum value of the two given values.
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
   Return the minimum value of the two given values.
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
   Return the absolute maximum value of the two given values.
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
   Return the absolute minimum value of the two given values.
*/
struct AbsMin : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return std::min(fabs(x),fabs(y));
  }
};

/**
   \brief
   Return the length between the two given values.
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
   Return the trace of the two given values.
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
   Perform the given operation upon all values in all data-points in the
   given Data object and return the final result.

   Calls DataArrayView::reductionOp
*/
template <class UnaryFunction>
inline
double
algorithm(DataExpanded& data,
          UnaryFunction operation)
{
  int i,j;
  int numDPPSample=data.getNumDPPSample();
  int numSamples=data.getNumSamples();
  int resultVectorLength=numDPPSample*numSamples;
  std::vector<double> resultVector(resultVectorLength);
  DataArrayView dataView=data.getPointDataView();
  // calculate the reduction operation value for each data point
  // storing the result for each data-point in successive entries
  // in resultVector
  {
#pragma omp for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
#pragma omp critical (reductionOp)
	resultVector[j*numSamples+i]=dataView.reductionOp(data.getPointOffset(i,j),operation);
      }
    }
  }
  // now calculate the reduction operation value across the results
  // for each data-point
  operation.resetResult();
  for (int l=0;l<resultVectorLength;l++) {
    operation(resultVector[l]);
  }
  return operation.getResult();
}

template <class UnaryFunction>
inline
double
algorithm(DataTagged& data,
          UnaryFunction operation)
{
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataArrayView& dataView=data.getPointDataView();
  std::vector<double> resultVector;
  int resultVectorLength;
  // perform the operation on each tagged value
  for (i=lookup.begin();i!=lookupEnd;i++) {
    resultVector.push_back(dataView.reductionOp(i->second,operation));
  }
  // perform the operation on the default value
  resultVector.push_back(data.getDefaultValue().reductionOp(operation));
  // now calculate the reduction operation value across the results
  // for each tagged value
  resultVectorLength=resultVector.size();
  operation.resetResult();
  for (int l=0;l<resultVectorLength;l++) {
    operation(resultVector[l]);
  }
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
   Perform the given data-point reduction operation on all data-points
   in data, storing results in corresponding data-points of result.

   Objects data and result must be of the same type, and have the same number
   of data points, but where data has data points of rank n, result must have
   data points of rank 0.

   Calls DataArrayView::reductionOp
*/
template <class UnaryFunction>
inline
void
dp_algorithm(DataExpanded& data,
             DataExpanded& result,
             UnaryFunction operation)
{
  int i,j;
  int numDPPSample=data.getNumDPPSample();
  int numSamples=data.getNumSamples();
  DataArrayView dataView=data.getPointDataView();
  DataArrayView resultView=result.getPointDataView();
  // perform the operation on each data-point and assign
  // this to the corresponding element in result
  {
#pragma omp for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
#pragma omp critical (reductionOp)
        resultView.getData(data.getPointOffset(i,j)) =
          dataView.reductionOp(data.getPointOffset(i,j),operation);
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
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataArrayView dataView=data.getPointDataView();
  DataArrayView resultView=result.getPointDataView();
  // perform the operation on each tagged data value
  // and assign this to the corresponding element in result
  for (i=lookup.begin();i!=lookupEnd;i++) {
    resultView.getData(i->second) =
      dataView.reductionOp(i->second,operation);
  }
  // perform the operation on the default data value
  // and assign this to the default element in result
  resultView.getData(0) =
    data.getDefaultValue().reductionOp(operation);
}

template <class UnaryFunction>
inline
void
dp_algorithm(DataConstant& data,
             DataConstant& result,
             UnaryFunction operation)
{
  // perform the operation on the data value
  // and assign this to the element in result
  result.getPointDataView().getData(0) =
    data.getPointDataView().reductionOp(operation);
}

} // end of namespace
#endif
