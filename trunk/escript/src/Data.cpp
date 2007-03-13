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
  m_protected=false;
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
  m_protected=false;
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
  m_protected=false;
#if defined DOPROF
  // create entry in global profiling table for this object
  profData = dataProfTable.newData();
#endif
}

Data::Data(const Data& inData)
{
  m_data=inData.m_data;
  m_protected=inData.isProtected();
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
  m_protected=false;
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
  m_protected=false;
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
  m_protected=false;
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
  m_protected=false;
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
  m_protected=false;
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
  m_protected=false;
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
  m_protected=false;
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
Data::setProtection() 
{ 
   m_protected=true;
}

bool
Data::isProtected() const 
{ 
   return m_protected;
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

Data
Data::oneOver() const
{
#if defined DOPROF
  profData->where++;
#endif
  return escript::unaryOp(*this,bind1st(divides<double>(),1.));
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
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
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
          numArray[make_tuple(dataPoint,i)]=dataPointView(i);
        }
      }
      if (dataPointRank==2) {
        for (int i=0; i<dataPointShape[0]; i++) {
          for (int j=0; j<dataPointShape[1]; j++) {
            numArray[make_tuple(dataPoint,i,j)] = dataPointView(i,j);
          }
        }
      }
      if (dataPointRank==3) {
        for (int i=0; i<dataPointShape[0]; i++) {
          for (int j=0; j<dataPointShape[1]; j++) {
            for (int k=0; k<dataPointShape[2]; k++) {
              numArray[make_tuple(dataPoint,i,j,k)]=dataPointView(i,j,k);
            }
          }
        }
      }
      if (dataPointRank==4) {
        for (int i=0; i<dataPointShape[0]; i++) {
          for (int j=0; j<dataPointShape[1]; j++) {
            for (int k=0; k<dataPointShape[2]; k++) {
              for (int l=0; l<dataPointShape[3]; l++) {
                numArray[make_tuple(dataPoint,i,j,k,l)]=dataPointView(i,j,k,l);
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
Data:: getValueOfDataPoint(int dataPointNo) 
{
  size_t length=0;
  int i, j, k, l;
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

  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       //
       // Check a valid sample number has been supplied
       if ((sampleNo >= getNumSamples()) || (sampleNo < 0 )) {
           throw DataException("Error - Data::convertToNumArray: invalid sampleNo.");
       }
              
       //
       // Check a valid data point number has been supplied
       if ((dataPointNoInSample >= getNumDataPointsPerSample()) || (dataPointNoInSample < 0)) {
           throw DataException("Error - Data::convertToNumArray: invalid dataPointNoInSample.");
       }
       // TODO: global error handling
       // create a view of the data if it is stored locally
       DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNoInSample);
		
       switch( dataPointRank ){
			case 0 :
				numArray[0] = dataPointView();
				break;
			case 1 :		
				for( i=0; i<dataPointShape[0]; i++ )
					numArray[i]=dataPointView(i);
				break;
			case 2 :		
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++)
						numArray[make_tuple(i,j)]=dataPointView(i,j);
				break;
			case 3 :		
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++)
							numArray[make_tuple(i,j,k)]=dataPointView(i,j,k);
				break;
			case 4 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++ )
							for( l=0; l<dataPointShape[3]; l++)
								numArray[make_tuple(i,j,k,l)]=dataPointView(i,j,k,l);
				break;
	}
  }
  //
  // return the array
  return numArray;

}
void
Data::setValueOfDataPointToArray(int dataPointNo, const boost::python::numeric::array num_array)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // check rank
  if (num_array.getrank()<getDataPointRank()) 
      throw DataException("Rank of numarray does not match Data object rank");

  //
  // check shape of num_array
  for (int i=0; i<getDataPointRank(); i++) {
    if (extract<int>(num_array.getshape()[i])!=getDataPointShape()[i])
       throw DataException("Shape of numarray does not match Data object rank");
  }
  //
  // make sure data is expanded:
  if (!isExpanded()) {
    expand();
  }
  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       m_data->copyToDataPoint(sampleNo, dataPointNoInSample,num_array);
  } else {
       m_data->copyToDataPoint(-1, 0,num_array);
  }
}

void
Data::setValueOfDataPoint(int dataPointNo, const double value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // make sure data is expanded:
  if (!isExpanded()) {
    expand();
  }
  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       m_data->copyToDataPoint(sampleNo, dataPointNoInSample,value);
  } else {
       m_data->copyToDataPoint(-1, 0,value);
  }
}

const 
boost::python::numeric::array
Data::getValueOfGlobalDataPoint(int procNo, int dataPointNo) 
{
  size_t length=0;
  int i, j, k, l, pos;
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

  // added for the MPI communication
  length=1;
  for( i=0; i<arrayRank; i++ ) length *= arrayShape[i];
  double *tmpData = new double[length];

  //
  // load the values for the data point into the numeric array.

	// updated for the MPI case
	if( get_MPIRank()==procNo ){
             if (getNumDataPointsPerSample()>0) {
                int sampleNo = dataPointNo/getNumDataPointsPerSample();
                int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
                //
                // Check a valid sample number has been supplied
                if ((sampleNo >= getNumSamples()) || (sampleNo < 0 )) {
                  throw DataException("Error - Data::convertToNumArray: invalid sampleNo.");
                }
              
                //
                // Check a valid data point number has been supplied
                if ((dataPointNoInSample >= getNumDataPointsPerSample()) || (dataPointNoInSample < 0)) {
                  throw DataException("Error - Data::convertToNumArray: invalid dataPointNoInSample.");
                }
                // TODO: global error handling
		// create a view of the data if it is stored locally
		DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNoInSample);
		
		// pack the data from the view into tmpData for MPI communication
		pos=0;
		switch( dataPointRank ){
			case 0 :
				tmpData[0] = dataPointView();
				break;
			case 1 :		
				for( i=0; i<dataPointShape[0]; i++ )
					tmpData[i]=dataPointView(i);
				break;
			case 2 :		
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++, pos++ )
						tmpData[pos]=dataPointView(i,j);
				break;
			case 3 :		
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++, pos++ )
							tmpData[pos]=dataPointView(i,j,k);
				break;
			case 4 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++ )
							for( l=0; l<dataPointShape[3]; l++, pos++ )
								tmpData[pos]=dataPointView(i,j,k,l);
				break;
		}
            }
	}
        #ifdef PASO_MPI	
        // broadcast the data to all other processes
	MPI_Bcast( tmpData, length, MPI_DOUBLE, procNo, get_MPIComm() );
        #endif

	// unpack the data
	switch( dataPointRank ){
		case 0 :
			numArray[0]=tmpData[0];
			break;
		case 1 :		
			for( i=0; i<dataPointShape[0]; i++ )
				numArray[i]=tmpData[i];
			break;
		case 2 :		
			for( i=0; i<dataPointShape[0]; i++ )
				for( j=0; j<dataPointShape[1]; j++ )
				   numArray[make_tuple(i,j)]=tmpData[i+j*dataPointShape[0]];
			break;
		case 3 :		
			for( i=0; i<dataPointShape[0]; i++ )
				for( j=0; j<dataPointShape[1]; j++ )
					for( k=0; k<dataPointShape[2]; k++ )
						numArray[make_tuple(i,j,k)]=tmpData[i+dataPointShape[0]*(j*+k*dataPointShape[1])];
			break;
		case 4 :
			for( i=0; i<dataPointShape[0]; i++ )
				for( j=0; j<dataPointShape[1]; j++ )
					for( k=0; k<dataPointShape[2]; k++ )
						for( l=0; l<dataPointShape[3]; l++ )
					        	numArray[make_tuple(i,j,k,l)]=tmpData[i+dataPointShape[0]*(j*+dataPointShape[1]*(k+l*dataPointShape[2]))];
			break;
	}

	delete [] tmpData;	
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
Data::erf() const
{
#if defined DOPROF
  profData->unary++;
#endif
#ifndef _WIN32
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::erf);
#else
  return Data();
#endif
}

Data
Data::asinh() const
{
#if defined DOPROF
  profData->unary++;
#endif
#ifndef _WIN32
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::asinh);
#else
  return Data();
#endif
}

Data
Data::acosh() const
{
#if defined DOPROF
  profData->unary++;
#endif
#ifndef _WIN32
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::acosh);
#else
  return Data();
#endif
}

Data
Data::atanh() const
{
#if defined DOPROF
  profData->unary++;
#endif
#ifndef _WIN32
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::atanh);
#else
  return Data();
#endif
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
  double localValue, globalValue;
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial absolute maximum value to zero

  AbsMax abs_max_func;
  localValue = algorithm(abs_max_func,0);
#ifdef PASO_MPI
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::Linf() const
{
  double localValue, globalValue;
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial absolute minimum value to max double
  AbsMin abs_min_func;
  localValue = algorithm(abs_min_func,numeric_limits<double>::max());

#ifdef PASO_MPI
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::sup() const
{
  double localValue, globalValue;
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  localValue = algorithm(fmax_func,numeric_limits<double>::max()*-1);
#ifdef PASO_MPI
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::inf() const
{
  double localValue, globalValue;
#if defined DOPROF
  profData->reduction1++;
#endif
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  localValue = algorithm(fmin_func,numeric_limits<double>::max());
#ifdef PASO_MPI
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

/* TODO */
/* global reduction */
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
Data::swapaxes(const int axis0, const int axis1) const
{
     int axis0_tmp,axis1_tmp;
     #if defined DOPROF
     profData->unary++;
     #endif
     DataArrayView::ShapeType s=getDataPointShape();
     DataArrayView::ShapeType ev_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=getDataPointRank();
     if (rank<2) {
        throw DataException("Error - Data::swapaxes argument must have at least rank 2.");
     }
     if (axis0<0 || axis0>rank-1) {
        throw DataException("Error - Data::swapaxes: axis0 must be between 0 and rank-1=" + rank-1);
     }
     if (axis1<0 || axis1>rank-1) {
         throw DataException("Error - Data::swapaxes: axis1 must be between 0 and rank-1=" + rank-1);
     }
     if (axis0 == axis1) {
         throw DataException("Error - Data::swapaxes: axis indices must be different.");
     }
     if (axis0 > axis1) {
         axis0_tmp=axis1;
         axis1_tmp=axis0;
     } else {
         axis0_tmp=axis0;
         axis1_tmp=axis1;
     }
     for (int i=0; i<rank; i++) {
       if (i == axis0_tmp) {
          ev_shape.push_back(s[axis1_tmp]); 
       } else if (i == axis1_tmp) {
          ev_shape.push_back(s[axis0_tmp]); 
       } else {
          ev_shape.push_back(s[i]); 
       }
     }
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->swapaxes(ev.m_data.get(), axis0_tmp, axis1_tmp);
     return ev;

}

Data
Data::symmetric() const
{
     #if defined DOPROF
        profData->unary++;
     #endif
     // check input
     DataArrayView::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        if(s[0] != s[1]) 
           throw DataException("Error - Data::symmetric can only be calculated for rank 2 object with equal first and second dimension.");
     }
     else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
           throw DataException("Error - Data::symmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
     }
     else {
        throw DataException("Error - Data::symmetric can only be calculated for rank 2 or 4 object.");
     }
     Data ev(0.,getDataPointShape(),getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->symmetric(ev.m_data.get());
     return ev;
}

Data
Data::nonsymmetric() const
{
     #if defined DOPROF
        profData->unary++;
     #endif
     // check input
     DataArrayView::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        if(s[0] != s[1]) 
           throw DataException("Error - Data::nonsymmetric can only be calculated for rank 2 object with equal first and second dimension.");
        DataArrayView::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->nonsymmetric(ev.m_data.get());
        return ev;
     }
     else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
           throw DataException("Error - Data::nonsymmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
        DataArrayView::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        ev_shape.push_back(s[2]);
        ev_shape.push_back(s[3]);
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->nonsymmetric(ev.m_data.get());
        return ev;
     }
     else {
        throw DataException("Error - Data::nonsymmetric can only be calculated for rank 2 or 4 object.");
     }
}

Data
Data::trace(int axis_offset) const
{
     #if defined DOPROF
        profData->unary++;
     #endif
     DataArrayView::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        DataArrayView::ShapeType ev_shape;
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     if (getDataPointRank()==3) {
        DataArrayView::ShapeType ev_shape;
        if (axis_offset==0) {
          int s2=s[2];
          ev_shape.push_back(s2);
        }
        else if (axis_offset==1) {
          int s0=s[0];
          ev_shape.push_back(s0);
        }
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     if (getDataPointRank()==4) {
        DataArrayView::ShapeType ev_shape;
        if (axis_offset==0) {
          ev_shape.push_back(s[2]);
          ev_shape.push_back(s[3]);
        }
        else if (axis_offset==1) {
          ev_shape.push_back(s[0]);
          ev_shape.push_back(s[3]);
        }
	else if (axis_offset==2) {
	  ev_shape.push_back(s[0]);
	  ev_shape.push_back(s[1]);
	}
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
	m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     else {
        throw DataException("Error - Data::trace can only be calculated for rank 2, 3 or 4 object.");
     }
}

Data
Data::transpose(int axis_offset) const
{
     #if defined DOPROF
     profData->unary++;
     #endif
     DataArrayView::ShapeType s=getDataPointShape();
     DataArrayView::ShapeType ev_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=getDataPointRank();
     if (axis_offset<0 || axis_offset>rank) {
        throw DataException("Error - Data::transpose must have 0 <= axis_offset <= rank=" + rank);
     }
     for (int i=0; i<rank; i++) {
       int index = (axis_offset+i)%rank;
       ev_shape.push_back(s[index]); // Append to new shape
     }
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->transpose(ev.m_data.get(), axis_offset);
     return ev;
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
Data::minGlobalDataPoint() const
{
  // NB: calc_minGlobalDataPoint( had to be split off from minGlobalDataPoint( as boost::make_tuple causes an
  // abort (for unknown reasons) if there are openmp directives with it in the
  // surrounding function

  int DataPointNo;
  int ProcNo;
  calc_minGlobalDataPoint(ProcNo,DataPointNo);
  return make_tuple(ProcNo,DataPointNo);
}

void
Data::calc_minGlobalDataPoint(int& ProcNo,
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

#ifdef PASO_MPI
	// determine the processor on which the minimum occurs
	next = temp.getDataPoint(lowi,lowj)();
	int lowProc = 0;
	double *globalMins = new double[get_MPISize()+1];
	int error = MPI_Gather ( &next, 1, MPI_DOUBLE, globalMins, 1, MPI_DOUBLE, 0, get_MPIComm() );
	
	if( get_MPIRank()==0 ){
		next = globalMins[lowProc];
		for( i=1; i<get_MPISize(); i++ )
			if( next>globalMins[i] ){
				lowProc = i;
				next = globalMins[i];
			}
	}
	MPI_Bcast( &lowProc, 1, MPI_DOUBLE, 0, get_MPIComm() );

	delete [] globalMins;
	ProcNo = lowProc;
#else
	ProcNo = 0;
#endif
  DataPointNo = lowj + lowi * numDPPSample;
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
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator+=(const boost::python::object& right)
{
  Data tmp(right,getFunctionSpace(),false);
  binaryOp(tmp,plus<double>());
  return (*this);
}

Data&
Data::operator-=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator-=(const boost::python::object& right)
{
  Data tmp(right,getFunctionSpace(),false);
  binaryOp(tmp,minus<double>());
  return (*this);
}

Data&
Data::operator*=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator*=(const boost::python::object& right)
{
  Data tmp(right,getFunctionSpace(),false);
  binaryOp(tmp,multiplies<double>());
  return (*this);
}

Data&
Data::operator/=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  binaryOp(right,divides<double>());
  return (*this);
}

Data&
Data::operator/=(const boost::python::object& right)
{
  Data tmp(right,getFunctionSpace(),false);
  binaryOp(tmp,divides<double>());
  return (*this);
}

Data
Data::rpowO(const boost::python::object& left) const
{
  Data left_d(left,*this);
  return left_d.powD(*this);
}

Data
Data::powO(const boost::python::object& right) const
{
  Data tmp(right,getFunctionSpace(),false);
  return powD(tmp);
}

Data
Data::powD(const Data& right) const
{
  Data result;
  if (getDataPointRank()<right.getDataPointRank()) {
     result.copy(right); 
     result.binaryOp(*this,escript::rpow);
  } else {
     result.copy(*this);
     result.binaryOp(right,(Data::BinaryDFunPtr)::pow);
  }
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
  if (left.getDataPointRank()<right.getDataPointRank()) {
     result.copy(right);
     result+=left;
  } else {
     result.copy(left);
     result+=right;
  }
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
  if (left.getDataPointRank()<right.getDataPointRank()) {
     result=right.neg();
     result+=left;
  } else {
     result.copy(left);
     result-=right;
  }
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
  if (left.getDataPointRank()<right.getDataPointRank()) {
     result.copy(right);
     result*=left;
  } else {
     result.copy(left);
     result*=right;
  }
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
  if (left.getDataPointRank()<right.getDataPointRank()) {
     result=right.oneOver();
     result*=left;
  } else {
     result.copy(left);
     result/=right;
  }
  return result;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const boost::python::object& right)
{
  return left+Data(right,left.getFunctionSpace(),false);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const boost::python::object& right)
{
  return left-Data(right,left.getFunctionSpace(),false);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const boost::python::object& right)
{
  return left*Data(right,left.getFunctionSpace(),false);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const boost::python::object& right)
{
  return left/Data(right,left.getFunctionSpace(),false);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const boost::python::object& left, const Data& right)
{
  return Data(left,right.getFunctionSpace(),false)+right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const boost::python::object& left, const Data& right)
{
  return Data(left,right.getFunctionSpace(),false)-right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const boost::python::object& left, const Data& right)
{
  return Data(left,right.getFunctionSpace(),false)*right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const boost::python::object& left, const Data& right)
{
  return Data(left,right.getFunctionSpace(),false)/right;
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

/* TODO */
/* global reduction */
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

/* TODO */
/* global reduction */
Data
Data::getSlice(const DataArrayView::RegionType& region) const
{
#if defined DOPROF
  profData->slicing++;
#endif
  return Data(*this,region);
}

/* TODO */
/* global reduction */
void
Data::setItemO(const boost::python::object& key,
               const boost::python::object& value)
{
  Data tempData(value,getFunctionSpace());
  setItemD(key,tempData);
}

/* TODO */
/* global reduction */
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

/* TODO */
/* global reduction */
void
Data::setSlice(const Data& value,
               const DataArrayView::RegionType& region)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
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

/* TODO */
/* global reduction */
void
Data::setTaggedValue(int tagKey,
                     const boost::python::object& value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
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

/* TODO */
/* global reduction */
void
Data::setTaggedValueFromCPP(int tagKey,
                            const DataArrayView& value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
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

/* TODO */
/* global reduction */
int
Data::getTagNumber(int dpno)
{
  return m_data->getTagNumber(dpno);
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
  vector<int> referenceNumbers(noSamples);
  for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
    referenceNumbers[sampleNo] = getFunctionSpace().getReferenceIDOfSample(sampleNo);
  }
  vector<int> tagNumbers(noSamples);
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
  vector<int> referenceNumbers(noSamples);
  for (int sampleNo=0; sampleNo<noSamples; sampleNo++) {
    archiveFile.read(reinterpret_cast<char *>(&referenceNumbers[sampleNo]),sizeof(int));
  }
  vector<int> tagNumbers(noSamples);
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
    if (referenceNumbers[sampleNo] != fspace.getReferenceIDOfSample(sampleNo)) {
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

Data
escript::C_GeneralTensorProduct(Data& arg_0,
                     Data& arg_1,
                     int axis_offset,
                     int transpose)
{
  // General tensor product: res(SL x SR) = arg_0(SL x SM) * arg_1(SM x SR)
  // SM is the product of the last axis_offset entries in arg_0.getShape().

  #if defined DOPROF
    // profData->binary++;
  #endif

  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z, arg_1_Z;
  if (arg_0.getFunctionSpace()!=arg_1.getFunctionSpace()) {
    if (arg_0.probeInterpolation(arg_1.getFunctionSpace())) {
      arg_0_Z = arg_0.interpolate(arg_1.getFunctionSpace());
      arg_1_Z = Data(arg_1);
    }
    else if (arg_1.probeInterpolation(arg_0.getFunctionSpace())) {
      arg_1_Z=arg_1.interpolate(arg_0.getFunctionSpace());
      arg_0_Z =Data(arg_0);
    }
    else {
      throw DataException("Error - C_GeneralTensorProduct: arguments have incompatible function spaces.");
    }
  } else {
      arg_0_Z = Data(arg_0);
      arg_1_Z = Data(arg_1);
  }
  // Get rank and shape of inputs
  int rank0 = arg_0_Z.getDataPointRank();
  int rank1 = arg_1_Z.getDataPointRank();
  DataArrayView::ShapeType shape0 = arg_0_Z.getDataPointShape();
  DataArrayView::ShapeType shape1 = arg_1_Z.getDataPointShape();

  // Prepare for the loops of the product and verify compatibility of shapes
  int start0=0, start1=0;
  if (transpose == 0)		{}
  else if (transpose == 1)	{ start0 = axis_offset; }
  else if (transpose == 2)	{ start1 = rank1-axis_offset; }
  else				{ throw DataException("C_GeneralTensorProduct: Error - transpose should be 0, 1 or 2"); }

  // Adjust the shapes for transpose
  DataArrayView::ShapeType tmpShape0;
  DataArrayView::ShapeType tmpShape1;
  for (int i=0; i<rank0; i++)	{ tmpShape0.push_back( shape0[(i+start0)%rank0] ); }
  for (int i=0; i<rank1; i++)	{ tmpShape1.push_back( shape1[(i+start1)%rank1] ); }

#if 0
  // For debugging: show shape after transpose
  char tmp[100];
  std::string shapeStr;
  shapeStr = "(";
  for (int i=0; i<rank0; i++)	{ sprintf(tmp, "%d,", tmpShape0[i]); shapeStr += tmp; }
  shapeStr += ")";
  cout << "C_GeneralTensorProduct: Shape of arg0 is " << shapeStr << endl;
  shapeStr = "(";
  for (int i=0; i<rank1; i++)	{ sprintf(tmp, "%d,", tmpShape1[i]); shapeStr += tmp; }
  shapeStr += ")";
  cout << "C_GeneralTensorProduct: Shape of arg1 is " << shapeStr << endl;
#endif

  // Prepare for the loops of the product
  int SL=1, SM=1, SR=1;
  for (int i=0; i<rank0-axis_offset; i++)	{
    SL *= tmpShape0[i];
  }
  for (int i=rank0-axis_offset; i<rank0; i++)	{
    if (tmpShape0[i] != tmpShape1[i-(rank0-axis_offset)]) {
      throw DataException("C_GeneralTensorProduct: Error - incompatible shapes");
    }
    SM *= tmpShape0[i];
  }
  for (int i=axis_offset; i<rank1; i++)		{
    SR *= tmpShape1[i];
  }

  // Define the shape of the output
  DataArrayView::ShapeType shape2;
  for (int i=0; i<rank0-axis_offset; i++) { shape2.push_back(tmpShape0[i]); } // First part of arg_0_Z
  for (int i=axis_offset; i<rank1; i++)   { shape2.push_back(tmpShape1[i]); } // Last part of arg_1_Z

  // Declare output Data object
  Data res;

  if      (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());	// DataConstant output
    double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[0]);
    double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[0]);
    double *ptr_2 = &((res.getPointDataView().getData())[0]);
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
  }
  else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

    // Prepare the DataConstant input
    DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

    // Borrow DataTagged input from Data object
    DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP_1 Programming error - casting to DataTagged."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());	// DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare offset into DataConstant
    int offset_0 = tmp_0->getPointOffset(0,0);
    double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
    // Get the views
    DataArrayView view_1 = tmp_1->getDefaultValue();
    DataArrayView view_2 = tmp_2->getDefaultValue();
    // Get the pointers to the actual data
    double *ptr_1 = &((view_1.getData())[0]);
    double *ptr_2 = &((view_2.getData())[0]);
    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_1.begin();i!=lookup_1.end();i++) {
      tmp_2->addTaggedValue(i->first,tmp_2->getDefaultValue());
      DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
      DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
      double *ptr_1 = &view_1.getData(0);
      double *ptr_2 = &view_2.getData(0);
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_1,dataPointNo_1;
    int numSamples_1 = arg_1_Z.getNumSamples();
    int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
    int offset_0 = tmp_0->getPointOffset(0,0);
    #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
    for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
      for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
        int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
        double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
        double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
        double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP_0 Programming error - casting to DataTagged."); }

    // Prepare the DataConstant input
    DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_0_Z.getFunctionSpace());	// DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare offset into DataConstant
    int offset_1 = tmp_1->getPointOffset(0,0);
    double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
    // Get the views
    DataArrayView view_0 = tmp_0->getDefaultValue();
    DataArrayView view_2 = tmp_2->getDefaultValue();
    // Get the pointers to the actual data
    double *ptr_0 = &((view_0.getData())[0]);
    double *ptr_2 = &((view_2.getData())[0]);
    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
      tmp_2->addTaggedValue(i->first,tmp_2->getDefaultValue());
      DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
      DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
      double *ptr_0 = &view_0.getData(0);
      double *ptr_2 = &view_2.getData(0);
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Borrow DataTagged input from Data object
    DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());
    res.tag();	// DataTagged output
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Get the views
    DataArrayView view_0 = tmp_0->getDefaultValue();
    DataArrayView view_1 = tmp_1->getDefaultValue();
    DataArrayView view_2 = tmp_2->getDefaultValue();
    // Get the pointers to the actual data
    double *ptr_0 = &((view_0.getData())[0]);
    double *ptr_1 = &((view_1.getData())[0]);
    double *ptr_2 = &((view_2.getData())[0]);
    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Merge the tags
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
      tmp_2->addTaggedValue(i->first,tmp_2->getDefaultValue()); // use tmp_2 to get correct shape
    }
    for (i=lookup_1.begin();i!=lookup_1.end();i++) {
      tmp_2->addTaggedValue(i->first,tmp_2->getDefaultValue());
    }
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_2=tmp_2->getTagLookup();
    for (i=lookup_2.begin();i!=lookup_2.end();i++) {
      DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
      DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
      DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
      double *ptr_0 = &view_0.getData(0);
      double *ptr_1 = &view_1.getData(0);
      double *ptr_2 = &view_2.getData(0);
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
      double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
        double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {

    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    int offset_1 = tmp_1->getPointOffset(0,0);
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
        double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
        double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }


  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
      double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
        double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
        double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
        double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else {
    throw DataException("Error - C_GeneralTensorProduct: unknown combination of inputs");
  }

  return res;
}

DataAbstract*
Data::borrowData() const
{
  return m_data.get();
}

/* Member functions specific to the MPI implementation */

void
Data::print() 
{
  int i,j;
  
  printf( "Data is %dX%d\n", getNumSamples(), getNumDataPointsPerSample() );
  for( i=0; i<getNumSamples(); i++ )
  {
    printf( "[%6d]", i );
    for( j=0; j<getNumDataPointsPerSample(); j++ )
      printf( "\t%10.7g", (getSampleData(i))[j] );
    printf( "\n" );
  }
}

int
Data::get_MPISize() const
{
	int error, size;
#ifdef PASO_MPI
	error = MPI_Comm_size( get_MPIComm(), &size );
#else
	size = 1;
#endif
	return size;
}

int
Data::get_MPIRank() const
{
	int error, rank;
#ifdef PASO_MPI
	error = MPI_Comm_rank( get_MPIComm(), &rank );
#else
	rank = 0;
#endif
	return rank;
}

MPI_Comm
Data::get_MPIComm() const
{ 
#ifdef PASO_MPI
	return MPI_COMM_WORLD;
#else
	return -1;
#endif
}

