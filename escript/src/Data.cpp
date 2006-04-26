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

#include "Data.h"

#include "DataExpanded.h"
#include "DataConstant.h"
#include "DataTagged.h"
#include "DataEmpty.h"
#include "DataArray.h"
#include "DataArrayView.h"
#include "DataProf.h"
#include "FunctionSpaceFactory.h"
#include "AbstractContinuousDomain.h"
#include "UnaryFuncs.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <functional>
#include <math.h>

#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/long.hpp>

using namespace std;
using namespace boost::python;
using namespace boost;
using namespace escript;

#if defined DOPROF
//
// global table of profiling data for all Data objects
DataProf dataProfTable;
#endif

Data::Data()
{
  //
  // Default data is type DataEmpty
  DataAbstract* temp=new DataEmpty();
  shared_ptr<DataAbstract> temp_data(temp);
  m_data=temp_data;
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(double value,
           const tuple& shape,
           const FunctionSpace& what,
           bool expanded)
{
  DataArrayView::ShapeType dataPointShape;
  for (int i = 0; i < shape.attr("__len__")(); ++i) {
    dataPointShape.push_back(extract<const int>(shape[i]));
  }
  DataArray temp(dataPointShape,value);
  initialise(temp.getView(),what,expanded);
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(double value,
	   const DataArrayView::ShapeType& dataPointShape,
	   const FunctionSpace& what,
           bool expanded)
{
  DataArray temp(dataPointShape,value);
  pair<int,int> dataShape=what.getDataShape();
  initialise(temp.getView(),what,expanded);
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const Data& inData)
{
  m_data=inData.m_data;
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const Data& inData,
           const DataArrayView::RegionType& region)
{
  //
  // Create Data which is a slice of another Data
  DataAbstract* tmp = inData.m_data->getSlice(region);
  shared_ptr<DataAbstract> temp_data(tmp);
  m_data=temp_data;
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const Data& inData,
           const FunctionSpace& functionspace)
{
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
  if (inData.getFunctionSpace()==functionspace) {
    m_data=inData.m_data;
  } else {
    #if defined DOPROF
    profData->interpolate++;
    #endif
    Data tmp(0,inData.getPointDataView().getShape(),functionspace,true);
    // Note: Must use a reference or pointer to a derived object
    // in order to get polymorphic behaviour. Shouldn't really
    // be able to create an instance of AbstractDomain but that was done
    // as a boost:python work around which may no longer be required.
    const AbstractDomain& inDataDomain=inData.getDomain();
    if  (inDataDomain==functionspace.getDomain()) {
      inDataDomain.interpolateOnDomain(tmp,inData);
    } else {
      inDataDomain.interpolateACross(tmp,inData);
    }
    m_data=tmp.m_data;
  }
}

Data::Data(const DataTagged::TagListType& tagKeys,
           const DataTagged::ValueListType & values,
           const DataArrayView& defaultValue,
           const FunctionSpace& what,
           bool expanded)
{
  DataAbstract* temp=new DataTagged(tagKeys,values,defaultValue,what);
  shared_ptr<DataAbstract> temp_data(temp);
  m_data=temp_data;
  if (expanded) {
    expand();
  }
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const numeric::array& value,
	   const FunctionSpace& what,
           bool expanded)
{
  initialise(value,what,expanded);
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const DataArrayView& value,
	   const FunctionSpace& what,
           bool expanded)
{
  initialise(value,what,expanded);
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const object& value,
	   const FunctionSpace& what,
           bool expanded)
{
  numeric::array asNumArray(value);
  initialise(asNumArray,what,expanded);
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const object& value,
           const Data& other)
{
  //
  // Create DataConstant using the given value and all other parameters
  // copied from other. If value is a rank 0 object this Data
  // will assume the point data shape of other.
  DataArray temp(value);
  if (temp.getView().getRank()==0) {
    //
    // Create a DataArray with the scalar value for all elements
    DataArray temp2(other.getPointDataView().getShape(),temp.getView()());
    initialise(temp2.getView(),other.getFunctionSpace(),false);
  } else {
    //
    // Create a DataConstant with the same sample shape as other
    initialise(temp.getView(),other.getFunctionSpace(),false);
  }
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::~Data()
{

}

escriptDataC
Data::getDataC()
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

escriptDataC
Data::getDataC() const
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

const boost::python::tuple
Data::getShapeTuple() const
{
  const DataArrayView::ShapeType& shape=getDataPointShape();
  switch(getDataPointRank()) {
     case 0:
        return make_tuple();
     case 1:
        return make_tuple(long_(shape[0]));
     case 2:
        return make_tuple(long_(shape[0]),long_(shape[1]));
     case 3:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]));
     case 4:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]),long_(shape[3]));
     default:
        throw DataException("Error - illegal Data rank.");
  }
}

void
Data::copy(const Data& other)
{
  //
  // Perform a deep copy
  {
    DataExpanded* temp=dynamic_cast<DataExpanded*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataExpanded copy
      DataAbstract* newData=new DataExpanded(*temp);
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  {
    DataTagged* temp=dynamic_cast<DataTagged*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataTagged copy
      DataAbstract* newData=new DataTagged(*temp);
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  {
    DataConstant* temp=dynamic_cast<DataConstant*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataConstant copy
      DataAbstract* newData=new DataConstant(*temp);
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  {
    DataEmpty* temp=dynamic_cast<DataEmpty*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataEmpty copy
      DataAbstract* newData=new DataEmpty();
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  throw DataException("Error - Copy not implemented for this Data type.");
}

void
Data::copyWithMask(const Data& other,
                   const Data& mask)
{
  Data mask1;
  Data mask2;

  mask1 = mask.wherePositive();
  mask2.copy(mask1);

  mask1 *= other;
  mask2 *= *this;
  mask2 = *this - mask2;

  *this = mask1 + mask2;
}

bool
Data::isExpanded() const
{
  DataExpanded* temp=dynamic_cast<DataExpanded*>(m_data.get());
  return (temp!=0);
}

bool
Data::isTagged() const
{
  DataTagged* temp=dynamic_cast<DataTagged*>(m_data.get());
  return (temp!=0);
}

bool
Data::isEmpty() const
{
  DataEmpty* temp=dynamic_cast<DataEmpty*>(m_data.get());
  return (temp!=0);
}

bool
Data::isConstant() const
{
  DataConstant* temp=dynamic_cast<DataConstant*>(m_data.get());
  return (temp!=0);
}

void
Data::expand()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataConst);
    shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else if (isTagged()) {
    DataTagged* tempDataTag=dynamic_cast<DataTagged*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataTag);
    shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else if (isExpanded()) {
    //
    // do nothing
  } else if (isEmpty()) {
    throw DataException("Error - Expansion of DataEmpty not possible.");
  } else {
    throw DataException("Error - Expansion not implemented for this Data type.");
  }
}

void
Data::tag()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataTagged(*tempDataConst);
    shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else if (isTagged()) {
    // do nothing
  } else if (isExpanded()) {
    throw DataException("Error - Creating tag data from DataExpanded not possible.");
  } else if (isEmpty()) {
    throw DataException("Error - Creating tag data from DataEmpty not possible.");
  } else {
    throw DataException("Error - Tagging not implemented for this Data type.");
  }
}

void
Data::reshapeDataPoint(const DataArrayView::ShapeType& shape) 
{
  m_data->reshapeDataPoint(shape);
}

Data
Data::wherePositive() const
{
#if defined DOPROF
  profData->where++;
#endif
  return escript::unaryOp(*this,bind2nd(greater<double>(),0.0));
}

Data
Data::whereNegative() const
{
#if defined DOPROF
  profData->where++;
#endif
  return escript::unaryOp(*this,bind2nd(less<double>(),0.0));
}

Data
Data::whereNonNegative() const
{
#if defined DOPROF
  profData->where++;
#endif
  return escript::unaryOp(*this,bind2nd(greater_equal<double>(),0.0));
}

Data
Data::whereNonPositive() const
{
#if defined DOPROF
  profData->where++;
#endif
  return escript::unaryOp(*this,bind2nd(less_equal<double>(),0.0));
}

Data
Data::whereZero(double tol) const
{
#if defined DOPROF
  profData->where++;
#endif
  Data dataAbs=abs();
  return escript::unaryOp(dataAbs,bind2nd(less_equal<double>(),tol));
}

Data
Data::whereNonZero(double tol) const
{
#if defined DOPROF
  profData->where++;
#endif
  Data dataAbs=abs();
  return escript::unaryOp(dataAbs,bind2nd(greater<double>(),tol));
}

Data
Data::interpolate(const FunctionSpace& functionspace) const
{
#if defined DOPROF
  profData->interpolate++;
#endif
  return Data(*this,functionspace);
}

bool
Data::probeInterpolation(const FunctionSpace& functionspace) const
{
  if (getFunctionSpace()==functionspace) {
    return true;
  } else {
    const AbstractDomain& domain=getDomain();
    if  (domain==functionspace.getDomain()) {
      return domain.probeInterpolationOnDomain(getFunctionSpace().getTypeCode(),functionspace.getTypeCode());
    } else {
      return domain.probeInterpolationACross(getFunctionSpace().getTypeCode(),functionspace.getDomain(),functionspace.getTypeCode());
    }
  }
}

Data
Data::gradOn(const FunctionSpace& functionspace) const
{
#if defined DOPROF
  profData->grad++;
#endif
  if (functionspace.getDomain()!=getDomain())
    throw DataException("Error - gradient cannot be calculated on different domains.");
  DataArrayView::ShapeType grad_shape=getPointDataView().getShape();
  grad_shape.push_back(functionspace.getDim());
  Data out(0.0,grad_shape,functionspace,true);
  getDomain().setToGradient(out,*this);
  return out;
}

Data
Data::grad() const
{
  return gradOn(escript::function(getDomain()));
}

int
Data::getDataPointSize() const
{
  return getPointDataView().noValues();
}

DataArrayView::ValueType::size_type
Data::getLength() const
{
  return m_data->getLength();
}

const DataArrayView::ShapeType&
Data::getDataPointShape() const
{
  return getPointDataView().getShape();
}

void
Data::fillFromNumArray(const boost::python::numeric::array num_array) 
{
  //
  // check rank
  if (num_array.getrank()<getDataPointRank()) 
      throw DataException("Rank of numarray does not match Data object rank");

  //
  // check shape of num_array
  for (int i=0; i<getDataPointRank(); i++) {
    if (extract<int>(num_array.getshape()[i+1])!=getDataPointShape()[i])
       throw DataException("Shape of numarray does not match Data object rank");
  }

  //
  // make sure data is expanded:
  if (!isExpanded()) {
    expand();
  }

  //
  // and copy over
  m_data->copyAll(num_array);
}

const
boost::python::numeric::array
Data::convertToNumArray()
{
  //
  // determine the total number of data points
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDataPointsPerSample();
  int numDataPoints = numSamples * numDataPointsPerSample;

  //
  // determine the rank and shape of each data point
  int dataPointRank = getDataPointRank();
  DataArrayView::ShapeType dataPointShape = getDataPointShape();

  //
  // create the numeric array to be returned
  boost::python::numeric::array numArray(0.0);

  //
  // the rank of the returned numeric array will be the rank of
  // the data points, plus one. Where the rank of the array is n,
  // the last n-1 dimensions will be equal to the shape of the
  // data points, whilst the first dimension will be equal to the
  // total number of data points. Thus the array will consist of
  // a serial vector of the data points.
  int arrayRank = dataPointRank + 1;
  DataArrayView::ShapeType arrayShape;
  arrayShape.push_back(numDataPoints);
  for (int d=0; d<dataPointRank; d++) {
     arrayShape.push_back(dataPointShape[d]);
  }

  //
  // resize the numeric array to the shape just calculated
  if (arrayRank==1) {
    numArray.resize(arrayShape[0]);
  }
  if (arrayRank==2) {
    numArray.resize(arrayShape[0],arrayShape[1]);
  }
  if (arrayRank==3) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2]);
  }
  if (arrayRank==4) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3]);
  }
  if (arrayRank==5) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3],arrayShape[4]);
  }

  //
  // loop through each data point in turn, loading the values for that data point
  // into the numeric array.
  int dataPoint = 0;
  for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
      DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNo);
      if (dataPointRank==0) {
        numArray[dataPoint]=dataPointView();
      }
      if (dataPointRank==1) {
        for (int i=0; i<dataPointShape[0]; i++) {
          numArray[dataPoint][i]=dataPointView(i);
        }
      }
      if (dataPointRank==2) {
        for (int i=0; i<dataPointShape[0]; i++) {
          for (int j=0; j<dataPointShape[1]; j++) {
            numArray[dataPoint][i][j] = dataPointView(i,j);
          }
        }
      }
      if (dataPointRank==3) {
        for (int i=0; i<dataPointShape[0]; i++) {
          for (int j=0; j<dataPointShape[1]; j++) {
            for (int k=0; k<dataPointShape[2]; k++) {
              numArray[dataPoint][i][j][k]=dataPointView(i,j,k);
            }
          }
        }
      }
      if (dataPointRank==4) {
        for (int i=0; i<dataPointShape[0]; i++) {
          for (int j=0; j<dataPointShape[1]; j++) {
            for (int k=0; k<dataPointShape[2]; k++) {
              for (int l=0; l<dataPointShape[3]; l++) {
                numArray[dataPoint][i][j][k][l]=dataPointView(i,j,k,l);
              }
            }
          }
        }
      }
      dataPoint++;
    }
  }

  //
  // return the loaded array
  return numArray;
}

const
boost::python::numeric::array
Data::convertToNumArrayFromSampleNo(int sampleNo)
{
  //
  // Check a valid sample number has been supplied
  if (sampleNo >= getNumSamples()) {
    throw DataException("Error - Data::convertToNumArray: invalid sampleNo.");
  }

  //
  // determine the number of data points per sample
  int numDataPointsPerSample = getNumDataPointsPerSample();

  //
  // determine the rank and shape of each data point
  int dataPointRank = getDataPointRank();
  DataArrayView::ShapeType dataPointShape = getDataPointShape();

  //
  // create the numeric array to be returned
  boost::python::numeric::array numArray(0.0);

  //
  // the rank of the returned numeric array will be the rank of
  // the data points, plus one. Where the rank of the array is n,
  // the last n-1 dimensions will be equal to the shape of the
  // data points, whilst the first dimension will be equal to the
  // total number of data points. Thus the array will consist of
  // a serial vector of the data points.
  int arrayRank = dataPointRank + 1;
  DataArrayView::ShapeType arrayShape;
  arrayShape.push_back(numDataPointsPerSample);
  for (int d=0; d<dataPointRank; d++) {
     arrayShape.push_back(dataPointShape[d]);
  }

  //
  // resize the numeric array to the shape just calculated
  if (arrayRank==1) {
    numArray.resize(arrayShape[0]);
  }
  if (arrayRank==2) {
    numArray.resize(arrayShape[0],arrayShape[1]);
  }
  if (arrayRank==3) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2]);
  }
  if (arrayRank==4) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3]);
  }
  if (arrayRank==5) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3],arrayShape[4]);
  }

  //
  // loop through each data point in turn, loading the values for that data point
  // into the numeric array.
  for (int dataPoint = 0; dataPoint < numDataPointsPerSample; dataPoint++) {
    DataArrayView dataPointView = getDataPoint(sampleNo, dataPoint);
    if (dataPointRank==0) {
      numArray[dataPoint]=dataPointView();
    }
    if (dataPointRank==1) {
      for (int i=0; i<dataPointShape[0]; i++) {
        numArray[dataPoint][i]=dataPointView(i);
      }
    }
    if (dataPointRank==2) {
      for (int i=0; i<dataPointShape[0]; i++) {
        for (int j=0; j<dataPointShape[1]; j++) {
          numArray[dataPoint][i][j] = dataPointView(i,j);
        }
      }
    }
    if (dataPointRank==3) {
      for (int i=0; i<dataPointShape[0]; i++) {
        for (int j=0; j<dataPointShape[1]; j++) {
          for (int k=0; k<dataPointShape[2]; k++) {
            numArray[dataPoint][i][j][k]=dataPointView(i,j,k);
          }
        }
      }
    }
    if (dataPointRank==4) {
      for (int i=0; i<dataPointShape[0]; i++) {
        for (int j=0; j<dataPointShape[1]; j++) {
          for (int k=0; k<dataPointShape[2]; k++) {
            for (int l=0; l<dataPointShape[3]; l++) {
              numArray[dataPoint][i][j][k][l]=dataPointView(i,j,k,l);
            }
          }
        }
      }
    }
  }

  //
  // return the loaded array
  return numArray;
}

const
boost::python::numeric::array
Data::convertToNumArrayFromDPNo(int sampleNo,
                                int dataPointNo)
{
  //
  // Check a valid sample number has been supplied
  if (sampleNo >= getNumSamples()) {
    throw DataException("Error - Data::convertToNumArray: invalid sampleNo.");
  }

  //
  // Check a valid data point number has been supplied
  if (dataPointNo >= getNumDataPointsPerSample()) {
    throw DataException("Error - Data::convertToNumArray: invalid dataPointNo.");
  }

  //
  // determine the rank and shape of each data point
  int dataPointRank = getDataPointRank();
  DataArrayView::ShapeType dataPointShape = getDataPointShape();

  //
  // create the numeric array to be returned
  boost::python::numeric::array numArray(0.0);

  //
  // the shape of the returned numeric array will be the same
  // as that of the data point
  int arrayRank = dataPointRank;
  DataArrayView::ShapeType arrayShape = dataPointShape;

  //
  // resize the numeric array to the shape just calculated
  if (arrayRank==0) {
    numArray.resize(1);
  }
  if (arrayRank==1) {
    numArray.resize(arrayShape[0]);
  }
  if (arrayRank==2) {
    numArray.resize(arrayShape[0],arrayShape[1]);
  }
  if (arrayRank==3) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2]);
  }
  if (arrayRank==4) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3]);
  }

  //
  // load the values for the data point into the numeric array.
  DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNo);
  if (dataPointRank==0) {
    numArray[0]=dataPointView();
  }
  if (dataPointRank==1) {
    for (int i=0; i<dataPointShape[0]; i++) {
      numArray[i]=dataPointView(i);
    }
  }
  if (dataPointRank==2) {
    for (int i=0; i<dataPointShape[0]; i++) {
      for (int j=0; j<dataPointShape[1]; j++) {
        numArray[i][j] = dataPointView(i,j);
      }
    }
  }
  if (dataPointRank==3) {
    for (int i=0; i<dataPointShape[0]; i++) {
      for (int j=0; j<dataPointShape[1]; j++) {
        for (int k=0; k<dataPointShape[2]; k++) {
          numArray[i][j][k]=dataPointView(i,j,k);
        }
      }
    }
  }
  if (dataPointRank==4) {
    for (int i=0; i<dataPointShape[0]; i++) {
      for (int j=0; j<dataPointShape[1]; j++) {
        for (int k=0; k<dataPointShape[2]; k++) {
          for (int l=0; l<dataPointShape[3]; l++) {
            numArray[i][j][k][l]=dataPointView(i,j,k,l);
          }
        }
      }
    }
  }

  //
  // return the loaded array
  return numArray;
}

boost::python::numeric::array
Data::integrate() const
{
  int index;
  int rank = getDataPointRank();
  DataArrayView::ShapeType shape = getDataPointShape();

#if defined DOPROF
  profData->integrate++;
#endif

  //
  // calculate the integral values
  vector<double> integrals(getDataPointSize());
  AbstractContinuousDomain::asAbstractContinuousDomain(getDomain()).setToIntegrals(integrals,*this);

  //
  // create the numeric array to be returned
  // and load the array with the integral values
  boost::python::numeric::array bp_array(1.0);
  if (rank==0) {
    bp_array.resize(1);
    index = 0;
    bp_array[0] = integrals[index];
  }
  if (rank==1) {
    bp_array.resize(shape[0]);
    for (int i=0; i<shape[0]; i++) {
      index = i;
      bp_array[i] = integrals[index];
    }
  }
  if (rank==2) {
       bp_array.resize(shape[0],shape[1]);
       for (int i=0; i<shape[0]; i++) {
         for (int j=0; j<shape[1]; j++) {
           index = i + shape[0] * j;
           bp_array[make_tuple(i,j)] = integrals[index];
         }
       }
  }
  if (rank==3) {
    bp_array.resize(shape[0],shape[1],shape[2]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        for (int k=0; k<shape[2]; k++) {
          index = i + shape[0] * ( j + shape[1] * k );
          bp_array[make_tuple(i,j,k)] = integrals[index];
        }
      }
    }
  }
  if (rank==4) {
    bp_array.resize(shape[0],shape[1],shape[2],shape[3]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        for (int k=0; k<shape[2]; k++) {
          for (int l=0; l<shape[3]; l++) {
            index = i + shape[0] * ( j + shape[1] * ( k + shape[2] * l ) );
            bp_array[make_tuple(i,j,k,l)] = integrals[index];
          }
        }
      }
    }
  }

  //
  // return the loaded array
  return bp_array;
}

Data
Data::sin() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::sin);
}

Data
Data::cos() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::cos);
}

Data
Data::tan() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::tan);
}

Data
Data::asin() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::asin);
}

Data
Data::acos() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::acos);
}

Data
Data::atan() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::atan);
}

Data
Data::sinh() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::sinh);
}

Data
Data::cosh() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::cosh);
}

Data
Data::tanh() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::tanh);
}

Data
Data::asinh() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::asinh);
}

Data
Data::acosh() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::acosh);
}

Data
Data::atanh() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::atanh);
}

Data
Data::log10() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::log10);
}

Data
Data::log() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::log);
}

Data
Data::sign() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,escript::fsign);
}

Data
Data::abs() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::fabs);
}

Data
Data::neg() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,negate<double>());
}

Data
Data::pos() const
{
#if defined DOPROF
  profData->unary++;
#endif
  Data result;
  // perform a deep copy
  result.copy(*this);
  return result;
}

Data
Data::exp() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::exp);
}

Data
Data::sqrt() const
{
#if defined DOPROF
  profData->unary++;
#endif
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::sqrt);
}

double
Data::Lsup() const
{
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial absolute maximum value to zero
  AbsMax abs_max_func;
  return algorithm(abs_max_func,0);
}

double
Data::Linf() const
{
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial absolute minimum value to max double
  AbsMin abs_min_func;
  return algorithm(abs_min_func,numeric_limits<double>::max());
}

double
Data::sup() const
{
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  return algorithm(fmax_func,numeric_limits<double>::max()*-1);
}

double
Data::inf() const
{
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  return algorithm(fmin_func,numeric_limits<double>::max());
}

Data
Data::maxval() const
{
#if defined DOPROF
  profData->reduction2++;
#endif
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  return dp_algorithm(fmax_func,numeric_limits<double>::max()*-1);
}

Data
Data::minval() const
{
#if defined DOPROF
  profData->reduction2++;
#endif
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  return dp_algorithm(fmin_func,numeric_limits<double>::max());
}

Data
Data::trace() const
{
#if defined DOPROF
  profData->reduction2++;
#endif
  Trace trace_func;
  return dp_algorithm(trace_func,0);
}

Data
Data::transpose(int axis) const
{
#if defined DOPROF
  profData->reduction2++;
#endif

  // not implemented
  throw DataException("Error - Data::transpose not implemented yet.");
  return Data();
}

Data
Data::eigenvalues() const
{
     #if defined DOPROF
        profData->unary++;
     #endif
     // check input
     DataArrayView::ShapeType s=getDataPointShape();
     if (getDataPointRank()!=2) 
        throw DataException("Error - Data::eigenvalues can only be calculated for rank 2 object.");
     if(s[0] != s[1]) 
        throw DataException("Error - Data::eigenvalues can only be calculated for object with equal first and second dimension.");
     // create return
     DataArrayView::ShapeType ev_shape(1,s[0]);
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->eigenvalues(ev.m_data.get());
     return ev;
}

const boost::python::tuple
Data::eigenvalues_and_eigenvectors(const double tol) const
{
     #if defined DOPROF
        profData->unary++;
     #endif
     DataArrayView::ShapeType s=getDataPointShape();
     if (getDataPointRank()!=2) 
        throw DataException("Error - Data::eigenvalues and eigenvectors can only be calculated for rank 2 object.");
     if(s[0] != s[1]) 
        throw DataException("Error - Data::eigenvalues and eigenvectors can only be calculated for object with equal first and second dimension.");
     // create return
     DataArrayView::ShapeType ev_shape(1,s[0]);
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     DataArrayView::ShapeType V_shape(2,s[0]);
     Data V(0.,V_shape,getFunctionSpace());
     V.typeMatchRight(*this);
     m_data->eigenvalues_and_eigenvectors(ev.m_data.get(),V.m_data.get(),tol);
     return make_tuple(boost::python::object(ev),boost::python::object(V));
}

const boost::python::tuple
Data::mindp() const
{
  // NB: calc_mindp had to be split off from mindp as boost::make_tuple causes an
  // abort (for unknown reasons) if there are openmp directives with it in the
  // surrounding function

  int SampleNo;
  int DataPointNo;

  calc_mindp(SampleNo,DataPointNo);

  return make_tuple(SampleNo,DataPointNo);
}

void
Data::calc_mindp(int& SampleNo,
                 int& DataPointNo) const
{
  int i,j;
  int lowi=0,lowj=0;
  double min=numeric_limits<double>::max();

  Data temp=minval();

  int numSamples=temp.getNumSamples();
  int numDPPSample=temp.getNumDataPointsPerSample();

  double next,local_min;
  int local_lowi,local_lowj;

  #pragma omp parallel private(next,local_min,local_lowi,local_lowj)
  {
    local_min=min;
    #pragma omp for private(i,j) schedule(static)
    for (i=0; i<numSamples; i++) {
      for (j=0; j<numDPPSample; j++) {
        next=temp.getDataPoint(i,j)();
        if (next<local_min) {
          local_min=next;
          local_lowi=i;
          local_lowj=j;
        }
      }
    }
    #pragma omp critical
    if (local_min<min) {
      min=local_min;
      lowi=local_lowi;
      lowj=local_lowj;
    }
  }

  SampleNo = lowi;
  DataPointNo = lowj;
}

void
Data::saveDX(std::string fileName) const
{
  boost::python::dict args;
  args["data"]=boost::python::object(this);
  getDomain().saveDX(fileName,args);
  return;
}

void
Data::saveVTK(std::string fileName) const
{
  boost::python::dict args;
  args["data"]=boost::python::object(this);
  getDomain().saveVTK(fileName,args);
  return;
}

Data&
Data::operator+=(const Data& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator+=(const boost::python::object& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator-=(const Data& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator-=(const boost::python::object& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator*=(const Data& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator*=(const boost::python::object& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator/=(const Data& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,divides<double>());
  return (*this);
}

Data&
Data::operator/=(const boost::python::object& right)
{
#if defined DOPROF
  profData->binary++;
#endif
  binaryOp(right,divides<double>());
  return (*this);
}

Data
Data::rpowO(const boost::python::object& left) const
{
#if defined DOPROF
  profData->binary++;
#endif
  Data left_d(left,*this);
  return left_d.powD(*this);
}

Data
Data::powO(const boost::python::object& right) const
{
#if defined DOPROF
  profData->binary++;
#endif
  Data result;
  result.copy(*this);
  result.binaryOp(right,(Data::BinaryDFunPtr)::pow);
  return result;
}

Data
Data::powD(const Data& right) const
{
#if defined DOPROF
  profData->binary++;
#endif
  Data result;
  result.copy(*this);
  result.binaryOp(right,(Data::BinaryDFunPtr)::pow);
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result+=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result-=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result*=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result/=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result+=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result-=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result*=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result/=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result+=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result-=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result*=right;
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result/=right;
  return result;
}

//
//bool escript::operator==(const Data& left, const Data& right)
//{
//  /*
//  NB: this operator does very little at this point, and isn't to
//  be relied on. Requires further implementation.
//  */
//
//  bool ret;
//
//  if (left.isEmpty()) {
//    if(!right.isEmpty()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
//  }
//
//  if (left.isConstant()) {
//    if(!right.isConstant()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
// }
//
//  if (left.isTagged()) {
//   if(!right.isTagged()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
//  }
//
//  if (left.isExpanded()) {
//    if(!right.isExpanded()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
//  }
//
//  return ret;
//}

Data
Data::getItem(const boost::python::object& key) const 
{
  const DataArrayView& view=getPointDataView();

  DataArrayView::RegionType slice_region=view.getSliceRegion(key);

  if (slice_region.size()!=view.getRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }

  return getSlice(slice_region);
}

Data
Data::getSlice(const DataArrayView::RegionType& region) const
{
#if defined DOPROF
  profData->slicing++;
#endif
  return Data(*this,region);
}

void
Data::setItemO(const boost::python::object& key,
               const boost::python::object& value)
{
  Data tempData(value,getFunctionSpace());
  setItemD(key,tempData);
}

void
Data::setItemD(const boost::python::object& key,
               const Data& value)
{
  const DataArrayView& view=getPointDataView();

  DataArrayView::RegionType slice_region=view.getSliceRegion(key);
  if (slice_region.size()!=view.getRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }
  if (getFunctionSpace()!=value.getFunctionSpace()) {
     setSlice(Data(value,getFunctionSpace()),slice_region);
  } else {
     setSlice(value,slice_region);
  }
}

void
Data::setSlice(const Data& value,
               const DataArrayView::RegionType& region)
{
#if defined DOPROF
  profData->slicing++;
#endif
  Data tempValue(value);
  typeMatchLeft(tempValue);
  typeMatchRight(tempValue);
  m_data->setSlice(tempValue.m_data.get(),region);
}

void
Data::typeMatchLeft(Data& right) const
{
  if (isExpanded()){
    right.expand();
  } else if (isTagged()) {
    if (right.isConstant()) {
      right.tag();
    }
  }
}

void
Data::typeMatchRight(const Data& right)
{
  if (isTagged()) {
    if (right.isExpanded()) {
      expand();
    }
  } else if (isConstant()) {
    if (right.isExpanded()) {
      expand();
    } else if (right.isTagged()) {
      tag();
    }
  }
}

void
Data::setTaggedValue(int tagKey,
                     const boost::python::object& value)
{
  //
  // Ensure underlying data object is of type DataTagged
  tag();

  if (!isTagged()) {
    throw DataException("Error - DataTagged conversion failed!!");
  }

  //
  // Construct DataArray from boost::python::object input value
  DataArray valueDataArray(value);

  //
  // Call DataAbstract::setTaggedValue
  m_data->setTaggedValue(tagKey,valueDataArray.getView());
}

void
Data::setTaggedValueFromCPP(int tagKey,
                            const DataArrayView& value)
{
  //
  // Ensure underlying data object is of type DataTagged
  tag();

  if (!isTagged()) {
    throw DataException("Error - DataTagged conversion failed!!");
  }
                                                                                                               
  //
  // Call DataAbstract::setTaggedValue
  m_data->setTaggedValue(tagKey,value);
}

int
Data::getTagNumber(int dpno)
{
  return m_data->getTagNumber(dpno);
}

void
Data::setRefValue(int ref,
                  const boost::python::numeric::array& value)
{
  //
  // Construct DataArray from boost::python::object input value
  DataArray valueDataArray(value);

  //
  // Call DataAbstract::setRefValue
  m_data->setRefValue(ref,valueDataArray);
}

void
Data::getRefValue(int ref,
                  boost::python::numeric::array& value)
{
  //
  // Construct DataArray for boost::python::object return value
  DataArray valueDataArray(value);

  //
  // Load DataArray with values from data-points specified by ref
  m_data->getRefValue(ref,valueDataArray);

  //
  // Load values from valueDataArray into return numarray

  // extract the shape of the numarray
  int rank = value.getrank();
  DataArrayView::ShapeType shape;
  for (int i=0; i < rank; i++) {
    shape.push_back(extract<int>(value.getshape()[i]));
  }

  // and load the numarray with the data from the DataArray
  DataArrayView valueView = valueDataArray.getView();

  if (rank==0) {
      boost::python::numeric::array temp_numArray(valueView());
      value = temp_numArray;
  }
  if (rank==1) {
    for (int i=0; i < shape[0]; i++) {
      value[i] = valueView(i);
    }
  }
  if (rank==2) {
    for (int i=0; i < shape[0]; i++) {
      for (int j=0; j < shape[1]; j++) {
        value[i][j] = valueView(i,j);
      }
    }
  }
  if (rank==3) {
    for (int i=0; i < shape[0]; i++) {
      for (int j=0; j < shape[1]; j++) {
        for (int k=0; k < shape[2]; k++) {
          value[i][j][k] = valueView(i,j,k);
        }
      }
    }
  }
  if (rank==4) {
    for (int i=0; i < shape[0]; i++) {
      for (int j=0; j < shape[1]; j++) {
        for (int k=0; k < shape[2]; k++) {
          for (int l=0; l < shape[3]; l++) {
            value[i][j][k][l] = valueView(i,j,k,l);
          }
        }
      }
    }
  }

}

void
Data::archiveData(const std::string fileName)
{
  cout << "Archiving Data object to: " << fileName << endl;

  //
  // Determine type of this Data object
  int dataType = -1;

  if (isEmpty()) {
    dataType = 0;
    cout << "\tdataType: DataEmpty" << endl;
  }
  if (isConstant()) {
    dataType = 1;
    cout << "\tdataType: DataConstant" << endl;
  }
  if (isTagged()) {
    dataType = 2;
    cout << "\tdataType: DataTagged" << endl;
  }
  if (isExpanded()) {
    dataType = 3;
    cout << "\tdataType: DataExpanded" << endl;
  }

  if (dataType == -1) {
    throw DataException("archiveData Error: undefined dataType");
  }

  //
  // Collect data items common to all Data types
  int noSamples = getNumSamples();
  int noDPPSample = getNumDataPointsPerSample();
  int functionSpaceType = getFunctionSpace().getTypeCode();
  int dataPointRank = getDataPointRank();
  int dataPointSize = getDataPointSize();
  int dataLength = getLength();
  DataArrayView::ShapeType dataPointShape = getDataPointShape();
  int referenceNumbers[noSamples];
  for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
    referenceNumbers[sampleNo] = getFunctionSpace().getReferenceNoFromSampleNo(sampleNo);
  }
  int tagNumbers[noSamples];
  if (isTagged()) {
    for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
      tagNumbers[sampleNo] = getFunctionSpace().getTagFromSampleNo(sampleNo);
    }
  }

  cout << "\tnoSamples: " << noSamples << " noDPPSample: " << noDPPSample << endl;
  cout << "\tfunctionSpaceType: " << functionSpaceType << endl;
  cout << "\trank: " << dataPointRank << " size: " << dataPointSize << " length: " << dataLength << endl;

  //
  // Flatten Shape to an array of integers suitable for writing to file
  int flatShape[4] = {0,0,0,0};
  cout << "\tshape: < ";
  for (int dim=0; dim<dataPointRank; dim++) {
    flatShape[dim] = dataPointShape[dim];
    cout << dataPointShape[dim] << " ";
  }
  cout << ">" << endl;

  //
  // Open archive file
  ofstream archiveFile;
  archiveFile.open(fileName.data(), ios::out);

  if (!archiveFile.good()) {
    throw DataException("archiveData Error: problem opening archive file");
  }

  //
  // Write common data items to archive file
  archiveFile.write(reinterpret_cast<char *>(&dataType),sizeof(int));
  archiveFile.write(reinterpret_cast<char *>(&noSamples),sizeof(int));
  archiveFile.write(reinterpret_cast<char *>(&noDPPSample),sizeof(int));
  archiveFile.write(reinterpret_cast<char *>(&functionSpaceType),sizeof(int));
  archiveFile.write(reinterpret_cast<char *>(&dataPointRank),sizeof(int));
  archiveFile.write(reinterpret_cast<char *>(&dataPointSize),sizeof(int));
  archiveFile.write(reinterpret_cast<char *>(&dataLength),sizeof(int));
  for (int dim = 0; dim < 4; dim++) {
    archiveFile.write(reinterpret_cast<char *>(&flatShape[dim]),sizeof(int));
  }
  for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
    archiveFile.write(reinterpret_cast<char *>(&referenceNumbers[sampleNo]),sizeof(int));
  }
  if (isTagged()) {
    for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
      archiveFile.write(reinterpret_cast<char *>(&tagNumbers[sampleNo]),sizeof(int));
    }
  }

  if (!archiveFile.good()) {
    throw DataException("archiveData Error: problem writing to archive file");
  }

  //
  // Archive underlying data values for each Data type
  int noValues;
  switch (dataType) {
    case 0:
      // DataEmpty
      noValues = 0;
      archiveFile.write(reinterpret_cast<char *>(&noValues),sizeof(int));
      cout << "\tnoValues: " << noValues << endl;
      break;
    case 1:
      // DataConstant
      noValues = m_data->getLength();
      archiveFile.write(reinterpret_cast<char *>(&noValues),sizeof(int));
      cout << "\tnoValues: " << noValues << endl;
      if (m_data->archiveData(archiveFile,noValues)) {
        throw DataException("archiveData Error: problem writing data to archive file");
      }
      break;
    case 2:
      // DataTagged
      noValues = m_data->getLength();
      archiveFile.write(reinterpret_cast<char *>(&noValues),sizeof(int));
      cout << "\tnoValues: " << noValues << endl;
      if (m_data->archiveData(archiveFile,noValues)) {
        throw DataException("archiveData Error: problem writing data to archive file");
      }
      break;
    case 3:
      // DataExpanded
      noValues = m_data->getLength();
      archiveFile.write(reinterpret_cast<char *>(&noValues),sizeof(int));
      cout << "\tnoValues: " << noValues << endl;
      if (m_data->archiveData(archiveFile,noValues)) {
        throw DataException("archiveData Error: problem writing data to archive file");
      }
      break;
  }

  if (!archiveFile.good()) {
    throw DataException("archiveData Error: problem writing data to archive file");
  }

  //
  // Close archive file
  archiveFile.close();

  if (!archiveFile.good()) {
    throw DataException("archiveData Error: problem closing archive file");
  }

}

void
Data::extractData(const std::string fileName,
                  const FunctionSpace& fspace)
{
  //
  // Can only extract Data to an object which is initially DataEmpty
  if (!isEmpty()) {
    throw DataException("extractData Error: can only extract to DataEmpty object");
  }

  cout << "Extracting Data object from: " << fileName << endl;

  int dataType;
  int noSamples;
  int noDPPSample;
  int functionSpaceType;
  int dataPointRank;
  int dataPointSize;
  int dataLength;
  DataArrayView::ShapeType dataPointShape;
  int flatShape[4];

  //
  // Open the archive file
  ifstream archiveFile;
  archiveFile.open(fileName.data(), ios::in);

  if (!archiveFile.good()) {
    throw DataException("extractData Error: problem opening archive file");
  }

  //
  // Read common data items from archive file
  archiveFile.read(reinterpret_cast<char *>(&dataType),sizeof(int));
  archiveFile.read(reinterpret_cast<char *>(&noSamples),sizeof(int));
  archiveFile.read(reinterpret_cast<char *>(&noDPPSample),sizeof(int));
  archiveFile.read(reinterpret_cast<char *>(&functionSpaceType),sizeof(int));
  archiveFile.read(reinterpret_cast<char *>(&dataPointRank),sizeof(int));
  archiveFile.read(reinterpret_cast<char *>(&dataPointSize),sizeof(int));
  archiveFile.read(reinterpret_cast<char *>(&dataLength),sizeof(int));
  for (int dim = 0; dim < 4; dim++) {
    archiveFile.read(reinterpret_cast<char *>(&flatShape[dim]),sizeof(int));
    if (flatShape[dim]>0) {
      dataPointShape.push_back(flatShape[dim]);
    }
  }
  int referenceNumbers[noSamples];
  for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
    archiveFile.read(reinterpret_cast<char *>(&referenceNumbers[sampleNo]),sizeof(int));
  }
  int tagNumbers[noSamples];
  if (dataType==2) {
    for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
      archiveFile.read(reinterpret_cast<char *>(&tagNumbers[sampleNo]),sizeof(int));
    }
  }

  if (!archiveFile.good()) {
    throw DataException("extractData Error: problem reading from archive file");
  }

  //
  // Verify the values just read from the archive file
  switch (dataType) {
    case 0:
      cout << "\tdataType: DataEmpty" << endl;
      break;
    case 1:
      cout << "\tdataType: DataConstant" << endl;
      break;
    case 2:
      cout << "\tdataType: DataTagged" << endl;
      break;
    case 3:
      cout << "\tdataType: DataExpanded" << endl;
      break;
    default:
      throw DataException("extractData Error: undefined dataType read from archive file");
      break;
  }

  cout << "\tnoSamples: " << noSamples << " noDPPSample: " << noDPPSample << endl;
  cout << "\tfunctionSpaceType: " << functionSpaceType << endl;
  cout << "\trank: " << dataPointRank << " size: " << dataPointSize << " length: " << dataLength << endl;
  cout << "\tshape: < ";
  for (int dim = 0; dim < dataPointRank; dim++) {
    cout << dataPointShape[dim] << " ";
  }
  cout << ">" << endl;

  //
  // Verify that supplied FunctionSpace object is compatible with this Data object.
  if ( (fspace.getTypeCode()!=functionSpaceType) ||
       (fspace.getNumSamples()!=noSamples) ||
       (fspace.getNumDPPSample()!=noDPPSample)
     ) {
    throw DataException("extractData Error: incompatible FunctionSpace");
  }
  for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
    if (referenceNumbers[sampleNo] != fspace.getReferenceNoFromSampleNo(sampleNo)) {
      throw DataException("extractData Error: incompatible FunctionSpace");
    }
  }
  if (dataType==2) {
    for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
      if (tagNumbers[sampleNo] != fspace.getTagFromSampleNo(sampleNo)) {
        throw DataException("extractData Error: incompatible FunctionSpace");
      }
    }
  }

  //
  // Construct a DataVector to hold underlying data values
  DataVector dataVec(dataLength);

  //
  // Load this DataVector with the appropriate values
  int noValues;
  archiveFile.read(reinterpret_cast<char *>(&noValues),sizeof(int));
  cout << "\tnoValues: " << noValues << endl;
  switch (dataType) {
    case 0:
      // DataEmpty
      if (noValues != 0) {
        throw DataException("extractData Error: problem reading data from archive file");
      }
      break;
    case 1:
      // DataConstant
      if (dataVec.extractData(archiveFile,noValues)) {
        throw DataException("extractData Error: problem reading data from archive file");
      }
      break;
    case 2:
      // DataTagged
      if (dataVec.extractData(archiveFile,noValues)) {
        throw DataException("extractData Error: problem reading data from archive file");
      }
      break;
    case 3:
      // DataExpanded
      if (dataVec.extractData(archiveFile,noValues)) {
        throw DataException("extractData Error: problem reading data from archive file");
      }
      break;
  }

  if (!archiveFile.good()) {
    throw DataException("extractData Error: problem reading from archive file");
  }

  //
  // Close archive file
  archiveFile.close();

  if (!archiveFile.good()) {
    throw DataException("extractData Error: problem closing archive file");
  }

  //
  // Construct an appropriate Data object
  DataAbstract* tempData;
  switch (dataType) {
    case 0:
      // DataEmpty
      tempData=new DataEmpty();
      break;
    case 1:
      // DataConstant
      tempData=new DataConstant(fspace,dataPointShape,dataVec);
      break;
    case 2:
      // DataTagged
      tempData=new DataTagged(fspace,dataPointShape,tagNumbers,dataVec);
      break;
    case 3:
      // DataExpanded
      tempData=new DataExpanded(fspace,dataPointShape,dataVec);
      break;
  }
  shared_ptr<DataAbstract> temp_data(tempData);
  m_data=temp_data;
}

ostream& escript::operator<<(ostream& o, const Data& data)
{
  o << data.toString();
  return o;
}
