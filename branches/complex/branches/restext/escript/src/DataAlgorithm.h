
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_DataAlgorithm_20040714_H
#define escript_DataAlgorithm_20040714_H
#include "system_dep.h"

#include "DataExpanded.h"
#include "DataTagged.h"
#include "DataConstant.h"

#include "DataMaths.h"

#include <iostream>
#include <algorithm>
#include <list>

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
    DataAlgorithmAdapter(const DataAlgorithmAdapter& other):
      m_initialValue(other.m_initialValue),
      m_currentValue(other.m_initialValue),
      operation(other.operation)
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
   \brief Return 1 if abs(x)>y, otherwise return 0.
*/
struct AbsGT : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return fabs(x)>y;
  }
};

/**
   \brief Return 1 if abs(x)<=y, otherwise return 0.
*/
struct AbsLTE : public std::binary_function<double,double,double>
{
  inline double operator()(double x, double y) const
  {
    return fabs(x)<=y;
  }
};


/**
   \brief
   Perform the given operation upon all values in all data-points in the
   given Data object and return the final result.
*/
template <class BinaryFunction>
inline
double
algorithm(const DataExpanded& data,
          BinaryFunction operation,
	  double initial_value)
{
  int i,j;
  int numDPPSample=data.getNumDPPSample();
  int numSamples=data.getNumSamples();
  double global_current_value=initial_value;
  double local_current_value;
//  DataArrayView dataView=data.getPointDataView();
  const DataTypes::ValueType& vec=data.getVectorRO();
  const DataTypes::ShapeType& shape=data.getShape();
  // calculate the reduction operation value for each data point
  // reducing the result for each data-point into the current_value variables
  #pragma omp parallel private(local_current_value)
  {
      local_current_value=initial_value;
      #pragma omp for private(i,j) schedule(static)
      for (i=0;i<numSamples;i++) {
        for (j=0;j<numDPPSample;j++) {
/*          local_current_value=operation(local_current_value,dataView.reductionOp(data.getPointOffset(i,j),operation,initial_value));*/
          local_current_value=operation(local_current_value,DataMaths::reductionOp(vec,shape,data.getPointOffset(i,j),operation,initial_value));

        }
      }
      #pragma omp critical
      global_current_value=operation(global_current_value,local_current_value);
  }
  return global_current_value;
}

// It is important that the algorithm only be applied to tags which are actually in use.
template <class BinaryFunction>
inline
double
algorithm(DataTagged& data,
          BinaryFunction operation,
	  double initial_value)
{
  double current_value=initial_value;

  const DataTypes::ValueType& vec=data.getVectorRO();
  const DataTypes::ShapeType& shape=data.getShape();
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  const std::list<int> used=data.getFunctionSpace().getListOfTagsSTL();
  for (std::list<int>::const_iterator i=used.begin();i!=used.end();++i)
  {
     int tag=*i;
     if (tag==0)	// check for the default tag
     {
	current_value=operation(current_value,DataMaths::reductionOp(vec,shape,data.getDefaultOffset(),operation,initial_value));
     }
     else
     {
	DataTagged::DataMapType::const_iterator it=lookup.find(tag);
	if (it!=lookup.end())
	{
		current_value=operation(current_value,DataMaths::reductionOp(vec,shape,it->second,operation,initial_value));
	}
     }
  }
  return current_value;
}

template <class BinaryFunction>
inline
double
algorithm(DataConstant& data,
          BinaryFunction operation,
	  double initial_value)
{
  return DataMaths::reductionOp(data.getVectorRO(),data.getShape(),0,operation,initial_value);
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
template <class BinaryFunction>
inline
void
dp_algorithm(const DataExpanded& data,
             DataExpanded& result,
             BinaryFunction operation,
	     double initial_value)
{
  int i,j;
  int numSamples=data.getNumSamples();
  int numDPPSample=data.getNumDPPSample();
//  DataArrayView dataView=data.getPointDataView();
//  DataArrayView resultView=result.getPointDataView();
  const DataTypes::ValueType& dataVec=data.getVectorRO();
  const DataTypes::ShapeType& shape=data.getShape();
  DataTypes::ValueType& resultVec=result.getVectorRW();
  // perform the operation on each data-point and assign
  // this to the corresponding element in result
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numSamples;i++) {
    for (j=0;j<numDPPSample;j++) {
/*      resultView.getData(result.getPointOffset(i,j)) =
        dataView.reductionOp(data.getPointOffset(i,j),operation,initial_value);*/
      resultVec[result.getPointOffset(i,j)] =
        DataMaths::reductionOp(dataVec, shape, data.getPointOffset(i,j),operation,initial_value);

    }
  }
}

template <class BinaryFunction>
inline
void
dp_algorithm(const DataTagged& data,
             DataTagged& result,
             BinaryFunction operation,
	     double initial_value)
{
  // perform the operation on each tagged value in data
  // and assign this to the corresponding element in result
  const DataTypes::ShapeType& shape=data.getShape();
  const DataTypes::ValueType& vec=data.getVectorRO();
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  for (DataTagged::DataMapType::const_iterator i=lookup.begin(); i!=lookup.end(); i++) {
    result.getDataByTagRW(i->first,0) =
	DataMaths::reductionOp(vec,shape,data.getOffsetForTag(i->first),operation,initial_value);
  }
  // perform the operation on the default data value
  // and assign this to the default element in result
  result.getVectorRW()[result.getDefaultOffset()] = DataMaths::reductionOp(data.getVectorRO(),data.getShape(),data.getDefaultOffset(),operation,initial_value);
}

template <class BinaryFunction>
inline
void
dp_algorithm(DataConstant& data,
             DataConstant& result,
             BinaryFunction operation,
	     double initial_value)
{
  // perform the operation on the data value
  // and assign this to the element in result
  result.getVectorRW()[0] =
    DataMaths::reductionOp(data.getVectorRO(),data.getShape(),0,operation,initial_value);
}

} // end of namespace

#endif
