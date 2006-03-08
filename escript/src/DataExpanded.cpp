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

#include "DataExpanded.h"
#include "DataException.h"
#include "DataConstant.h"
#include "DataTagged.h"

#include <boost/python/extract.hpp>

using namespace std;
using namespace boost::python;
using namespace boost;

namespace escript {

DataExpanded::DataExpanded(const boost::python::numeric::array& value,
                           const FunctionSpace& what)
  : DataAbstract(what)
{
  DataArrayView::ShapeType tempShape;
  //
  // extract the shape of the python numarray
  for (int i=0; i<value.getrank(); i++) {
    tempShape.push_back(extract<int>(value.getshape()[i]));
  }
  //
  // initialise the data array for this object
  initialise(tempShape,what.getNumSamples(),what.getNumDPPSample());
  //
  // copy the given value to every data point
  copy(value);
}

DataExpanded::DataExpanded(const DataExpanded& other)
  : DataAbstract(other.getFunctionSpace()),
  m_data(other.m_data)
{
  //
  // create the view for the data
  DataArrayView temp(m_data.getData(),other.getPointDataView().getShape());
  setPointDataView(temp);
}

DataExpanded::DataExpanded(const DataConstant& other)
  : DataAbstract(other.getFunctionSpace())
{
  //
  // initialise the data array for this object
  initialise(other.getPointDataView().getShape(),other.getNumSamples(),other.getNumDPPSample());
  //
  // DataConstant only has one value, copy this to every data point
  copy(other.getPointDataView());
}

DataExpanded::DataExpanded(const DataTagged& other)
  : DataAbstract(other.getFunctionSpace())
{
  //
  // initialise the data array for this object
  initialise(other.getPointDataView().getShape(),other.getNumSamples(),other.getNumDPPSample());
  //
  // for each data point in this object, extract and copy the corresponding data
  // value from the given DataTagged object
  int i,j;
  DataArrayView::ValueType::size_type numRows=m_data.getNumRows();
  DataArrayView::ValueType::size_type numCols=m_data.getNumCols();
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numRows;i++) {
    for (j=0;j<numCols;j++) {
      try {
        getPointDataView().copy(getPointOffset(i,j),
                                other.getPointDataView(),
                                other.getPointOffset(i,j));
      }
      catch (std::exception& e) {
        cout << e.what() << endl;
      }
    }
  }
}

DataExpanded::DataExpanded(const DataExpanded& other,
                           const DataArrayView::RegionType& region)
  : DataAbstract(other.getFunctionSpace())
{
  //
  // get the shape of the slice
  DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
  //
  // initialise this Data object to the shape of the slice
  initialise(shape,other.getNumSamples(),other.getNumDPPSample());
  //
  // copy the data
  DataArrayView::RegionLoopRangeType region_loop_range=getSliceRegionLoopRange(region);
  DataArrayView::ValueType::size_type numRows=m_data.getNumRows();
  DataArrayView::ValueType::size_type numCols=m_data.getNumCols();
  int i,j;
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numRows;i++) {
    for (j=0;j<numCols;j++) {
      try {
        getPointDataView().copySlice(getPointOffset(i,j),
                                     other.getPointDataView(),
                                     other.getPointOffset(i,j),
                                     region_loop_range);
      }
      catch (std::exception& e) {
        cout << e.what() << endl;
      }
    }
  }
}

DataExpanded::DataExpanded(const DataArrayView& value,
                           const FunctionSpace& what)
  : DataAbstract(what)
{
  //
  // get the shape of the given data value
  DataArrayView::ShapeType tempShape=value.getShape();
  //
  // initialise this Data object to the shape of the given data value
  initialise(tempShape,what.getNumSamples(),what.getNumDPPSample());
  //
  // copy the given value to every data point
  copy(value);
}

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataArrayView::ShapeType &shape,
                           const DataArrayView::ValueType &data)
  : DataAbstract(what)
{
  //
  // create the view of the data
  initialise(shape,what.getNumSamples(),what.getNumDPPSample());
  //
  // copy the data in the correct format
  m_data.getData()=data;
}

DataExpanded::~DataExpanded()
{
}

void
DataExpanded::reshapeDataPoint(const DataArrayView::ShapeType& shape) 
{
  if (getPointDataView().getRank()!=0) {
    stringstream temp;
    temp << "Error - Can only reshape Data with data points of rank 0. "
         << "This Data has data points with rank: " 
         << getPointDataView().getRank();
    throw DataException(temp.str());
  }
  //
  // create the new DataBlocks2D data array, and a corresponding DataArrayView
  DataBlocks2D newData(getNumSamples(),getNumDPPSample(),DataArrayView::noValues(shape));
  DataArrayView newView(newData.getData(),shape);
  //
  // Copy the original data to every value for the new shape
  int i,j;
  int nRows=m_data.getNumRows();
  int nCols=m_data.getNumCols();
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<nRows;i++) {
    for (j=0;j<nCols;j++) {
      // NOTE: An exception may be thown from this call if 
      // DOASSERT is on which of course will play
      // havoc with the omp threads. Run single threaded
      // if using DOASSERT. 
      newView.copy(newData.index(i,j),m_data(i,j));
    }
  }
  // swap the new data array for the original
  m_data.Swap(newData);
  // set the corresponding DataArrayView
  DataArrayView temp(m_data.getData(),shape);
  setPointDataView(temp);
}

DataAbstract*
DataExpanded::getSlice(const DataArrayView::RegionType& region) const 
{
  return new DataExpanded(*this,region);
}

void
DataExpanded::setSlice(const DataAbstract* value,
                       const DataArrayView::RegionType& region) 
{
  const DataExpanded* tempDataExp=dynamic_cast<const DataExpanded*>(value);
  if (tempDataExp==0) {
    throw DataException("Programming error - casting to DataExpanded.");
  }
  //
  // get shape of slice
  DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
  DataArrayView::RegionLoopRangeType region_loop_range=getSliceRegionLoopRange(region);
  //
  // check shape
  if (getPointDataView().getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (tempDataExp->getPointDataView().getRank()>0 and !value->getPointDataView().checkShape(shape)) {
    throw DataException (value->getPointDataView().createShapeErrorMessage(
		"Error - Couldn't copy slice due to shape mismatch.",shape));
  }
  //
  // copy the data from the slice into this object
  DataArrayView::ValueType::size_type numRows=m_data.getNumRows();
  DataArrayView::ValueType::size_type numCols=m_data.getNumCols();
  int i, j;
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numRows;i++) {
    for (j=0;j<numCols;j++) {
      getPointDataView().copySliceFrom(getPointOffset(i,j),
                                       tempDataExp->getPointDataView(),
                                       tempDataExp->getPointOffset(i,j),
                                       region_loop_range);
    }
  }
}

void
DataExpanded::copy(const DataArrayView& value) 
{
  //
  // copy a single value to every data point in this object
  int nRows=m_data.getNumRows();
  int nCols=m_data.getNumCols();
  int i,j;
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<nRows;i++) {
    for (j=0;j<nCols;j++) {
      // NOTE: An exception may be thown from this call if 
      // DOASSERT is on which of course will play
      // havoc with the omp threads. Run single threaded
      // if using DOASSERT. 
      getPointDataView().copy(m_data.index(i,j),value);
    }
  }
}

void
DataExpanded::copy(const boost::python::numeric::array& value) 
{
  //
  // first convert the numarray into a DataArray object
  DataArray temp(value);
  //
  // check the input shape matches this shape
  if (!getPointDataView().checkShape(temp.getView().getShape())) {
    throw DataException(getPointDataView().createShapeErrorMessage(
                        "Error - (DataExpanded) Cannot copy due to shape mismatch.",
                        temp.getView().getShape()));
  }
  //
  // now copy over the data
  copy(temp.getView());
}

void
DataExpanded::initialise(const DataArrayView::ShapeType& shape,
                         int noSamples,
                         int noDataPointsPerSample)
{
  //
  // resize data array to the required size
  m_data.resize(noSamples,noDataPointsPerSample,DataArrayView::noValues(shape));
  //
  // create the data view of the data array
  DataArrayView temp(m_data.getData(),shape);
  setPointDataView(temp);
}

string
DataExpanded::toString() const
{
  stringstream temp;
  //
  // create a temporary view as the offset will be changed
  DataArrayView tempView(getPointDataView().getData(),getPointDataView().getShape(),getPointDataView().getOffset());
  for (int i=0;i<m_data.getNumRows();i++) {
    for (int j=0;j<m_data.getNumCols();j++) {
      tempView.setOffset(m_data.index(i,j));
      stringstream suffix;
      suffix << "(" << i << "," << j << ")";
      temp << tempView.toString(suffix.str());
      if (!(i==(m_data.getNumRows()-1) && j==(m_data.getNumCols()-1))) {
        temp << endl;
      }
    }
  }
  return temp.str();
}

DataArrayView::ValueType::size_type
DataExpanded::getPointOffset(int sampleNo,
                             int dataPointNo) const
{
  return m_data.index(sampleNo,dataPointNo);
}

DataArrayView
DataExpanded::getDataPoint(int sampleNo,
                           int dataPointNo)
{
  DataArrayView temp(m_data.getData(),getPointDataView().getShape(),m_data.index(sampleNo,dataPointNo));
  return temp;
}

DataArrayView::ValueType::size_type
DataExpanded::getLength() const
{
  return m_data.size();
}

void
DataExpanded::setRefValue(int ref,
                          const DataArray& value)
{
  //
  // Get the number of samples and data-points per sample.
  int numSamples = getNumSamples();
  int numDPPSample = getNumDPPSample();

  //
  // Determine the sample number which corresponds to this reference number.
  int sampleNo = -1;
  int tempRef = -1;
  for (int n=0; n<numSamples; n++) {
    tempRef = getFunctionSpace().getReferenceNoFromSampleNo(n);
    if (tempRef == ref) {
      sampleNo = n;
      break;
    }
  }
  if (sampleNo == -1) {
    throw DataException("DataExpanded::setRefValue error: invalid ref number supplied.");
  }

  for (int n=0; n<numDPPSample; n++) {
    //
    // Get *each* data-point in the sample in turn.
    DataArrayView pointView = getDataPoint(sampleNo, n);
    //
    // Assign the values in the DataArray to this data-point.
    pointView.copy(value.getView());
  }
}

void
DataExpanded::getRefValue(int ref,
                          DataArray& value)
{
  //
  // Get the number of samples and data-points per sample.
  int numSamples = getNumSamples();
  int numDPPSample = getNumDPPSample();

  //
  // Determine the sample number which corresponds to this reference number
  int sampleNo = -1;
  int tempRef = -1;
  for (int n=0; n<numSamples; n++) {
    tempRef = getFunctionSpace().getReferenceNoFromSampleNo(n);
    if (tempRef == ref) {
      sampleNo = n;
      break;
    }
  }
  if (sampleNo == -1) {
    throw DataException("DataExpanded::getRefValue error: invalid ref number supplied.");
  }

  //
  // Get the *first* data-point associated with this sample number.
  DataArrayView pointView = getDataPoint(sampleNo, 0);

  //
  // Load the values from this data-point into the DataArray
  value.getView().copy(pointView);
}

int
DataExpanded::archiveData(ofstream& archiveFile,
                          const DataArrayView::ValueType::size_type noValues) const
{
  return(m_data.archiveData(archiveFile, noValues));
}

int
DataExpanded::extractData(ifstream& archiveFile,
                          const DataArrayView::ValueType::size_type noValues)
{
  return(m_data.extractData(archiveFile, noValues));
}

void
DataExpanded::copyAll(const boost::python::numeric::array& value) {
  //
  // Get the number of samples and data-points per sample.
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int dataPointRank = getPointDataView().getRank();
  ShapeType dataPointShape = getPointDataView().getShape();
  //
  // check rank:
  if (value.getrank()!=dataPointRank+1)
       throw DataException("Rank of numarray does not match Data object rank");
  if (value.getshape()[0]!=numSamples*numDataPointsPerSample)
       throw DataException("leading dimension of numarray is too small");
  //
  int dataPoint = 0;
  for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
      DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNo);
      if (dataPointRank==0) {
         dataPointView()=extract<double>(value[dataPoint]);
      } else if (dataPointRank==1) {
         for (int i=0; i<dataPointShape[0]; i++) {
            dataPointView(i)=extract<double>(value[dataPoint][i]);
         }
      } else if (dataPointRank==2) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             dataPointView(i,j)=extract<double>(value[dataPoint][i][j]);
           }
         }
       } else if (dataPointRank==3) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
                 dataPointView(i,j,k)=extract<double>(value[dataPoint][i][j][k]);
             }
           }
         }
       } else if (dataPointRank==4) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
               for (int l=0; l<dataPointShape[3]; l++) {
                 dataPointView(i,j,k,l)=extract<double>(value[dataPoint][i][j][k][l]);
               }
             }
           }
         }
      }
      dataPoint++;
    }
  }
}
void
DataExpanded::eigenvalues(DataAbstract* ev)
{
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::eigenvalues: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::eigenvalues(thisView,getPointOffset(sampleNo,dataPointNo),
                                    evView,ev->getPointOffset(sampleNo,dataPointNo));
    }
  }
}
void
DataExpanded::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::eigenvalues_and_eigenvectors: casting to DataExpanded failed (propably a programming error).");
  }
  DataExpanded* temp_V=dynamic_cast<DataExpanded*>(V);
  if (temp_V==0) {
    throw DataException("Error - DataExpanded::eigenvalues_and_eigenvectors: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView& VView=V->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::eigenvalues_and_eigenvectors(thisView,getPointOffset(sampleNo,dataPointNo),
                                                     evView,ev->getPointOffset(sampleNo,dataPointNo),
                                                     VView,V->getPointOffset(sampleNo,dataPointNo),
                                                     tol);
    }
  }
}

}  // end of namespace
