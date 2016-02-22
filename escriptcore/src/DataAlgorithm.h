
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


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
#include <cmath>	// for max

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
template<typename T>
struct AbsMax 
{
  inline DataTypes::real_t operator()(T x, T y) const
  {
    return std::max(std::abs(x),std::abs(y));
  }
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef DataTypes::real_t result_type;
};

/**
   \brief
   Return the absolute minimum value of the two given values.
*/
template<typename T>
struct AbsMin 
{
  inline DataTypes::real_t operator()(T x, T y) const
  {
    return min(abs(x),abs(y));
  }
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef DataTypes::real_t result_type;
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
DataTypes::real_t
algorithm(const DataExpanded& data,
          BinaryFunction operation,
	  DataTypes::real_t initial_value)
{
  int i,j;
  int numDPPSample=data.getNumDPPSample();
  int numSamples=data.getNumSamples();
  DataTypes::real_t global_current_value=initial_value;
  DataTypes::real_t local_current_value;
//  DataArrayView dataView=data.getPointDataView();
  const auto& vec=data.getTypedVectorRO(typename BinaryFunction::first_argument_type(0));
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
DataTypes::real_t
algorithm(DataTagged& data,
          BinaryFunction operation,
	  DataTypes::real_t initial_value)
{
  DataTypes::real_t current_value=initial_value;

  const auto& vec=data.getTypedVectorRO(typename BinaryFunction::first_argument_type(0));
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
DataTypes::real_t
algorithm(DataConstant& data,
          BinaryFunction operation,
	  DataTypes::real_t initial_value)
{
  return DataMaths::reductionOp(data.getTypedVectorRO(typename BinaryFunction::first_argument_type(0)),data.getShape(),0,operation,initial_value);
}

/**
   \brief
   Perform the given data-point reduction operation on all data-points
   in data, storing results in corresponding data-points of result.

   Objects data and result must be of the same type, and have the same number
   of data points, but where data has data points of rank n, result must have
   data points of rank 0.
*/
template <class BinaryFunction>
inline
void
dp_algorithm(const DataExpanded& data,
             DataExpanded& result,
             BinaryFunction operation,
	     typename BinaryFunction::first_argument_type initial_value)
{
  int i,j;
  int numSamples=data.getNumSamples();
  int numDPPSample=data.getNumDPPSample();
//  DataArrayView dataView=data.getPointDataView();
//  DataArrayView resultView=result.getPointDataView();
  const auto& dataVec=data.getTypedVectorRO(initial_value);
  const DataTypes::ShapeType& shape=data.getShape();
  auto& resultVec=result.getTypedVectorRW(initial_value);
  // perform the operation on each data-point and assign
  // this to the corresponding element in result
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numSamples;i++) {
    for (j=0;j<numDPPSample;j++) {
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
	     typename BinaryFunction::first_argument_type initial_value)
{
  // perform the operation on each tagged value in data
  // and assign this to the corresponding element in result
  const DataTypes::ShapeType& shape=data.getShape();
  const auto& vec=data.getTypedVectorRO(initial_value);
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  for (DataTagged::DataMapType::const_iterator i=lookup.begin(); i!=lookup.end(); i++) {
    result.getDataByTagRW(i->first,0) =
	DataMaths::reductionOp(vec,shape,data.getOffsetForTag(i->first),operation,initial_value);
  }
  // perform the operation on the default data value
  // and assign this to the default element in result
  result.getTypedVectorRW(initial_value)[result.getDefaultOffset()] = DataMaths::reductionOp(data.getTypedVectorRO(initial_value),data.getShape(),data.getDefaultOffset(),operation,initial_value);
}

template <class BinaryFunction>
inline
void
dp_algorithm(DataConstant& data,
             DataConstant& result,
             BinaryFunction operation,
	     typename BinaryFunction::first_argument_type initial_value)
{
  // perform the operation on the data value
  // and assign this to the element in result
  result.getTypedVectorRW(initial_value)[0] =
    DataMaths::reductionOp(data.getTypedVectorRO(initial_value),data.getShape(),0,operation,initial_value);
}

} // end of namespace

#endif
