
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "DataExpanded.h"
#include "DataException.h"
#include "DataConstant.h"
#include "DataTagged.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

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
  if (tempDataExp->getPointDataView().getRank()>0 && !value->getPointDataView().checkShape(shape)) {
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
      getPointDataView().copy(getPointOffset(i,j),value);
    }
  }
}

void
DataExpanded::copy(const boost::python::numeric::array& value)
{

  // extract the shape of the numarray
  DataArrayView::ShapeType tempShape;
  for (int i=0; i < value.getrank(); i++) {
    tempShape.push_back(extract<int>(value.getshape()[i]));
  }

  // get the space for the data vector
  int len = DataArrayView::noValues(tempShape);
  DataVector temp_data(len, 0.0, len);
  DataArrayView temp_dataView(temp_data, tempShape);
  temp_dataView.copy(value);

  //
  // check the input shape matches this shape
  if (!getPointDataView().checkShape(temp_dataView.getShape())) {
    throw DataException(getPointDataView().createShapeErrorMessage(
                        "Error - (DataExpanded) Cannot copy due to shape mismatch.",
                        temp_dataView.getShape()));
  }
  //
  // now copy over the data
  copy(temp_dataView);
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
  FunctionSpace fs=getFunctionSpace();
  //
  // create a temporary view as the offset will be changed
  DataArrayView tempView(getPointDataView().getData(),getPointDataView().getShape(),getPointDataView().getOffset());
  for (int i=0;i<m_data.getNumRows();i++) {
    for (int j=0;j<m_data.getNumCols();j++) {
      tempView.setOffset(getPointOffset(i,j));
      stringstream suffix;
      suffix << "( id: " << i << ", ref: " << fs.getReferenceIDOfSample(i) << ", pnt: " << j << ")";
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
  DataArrayView temp(m_data.getData(),getPointDataView().getShape(),getPointOffset(sampleNo,dataPointNo));
  return temp;
}

DataArrayView::ValueType::size_type
DataExpanded::getLength() const
{
  return m_data.size();
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
DataExpanded::copyToDataPoint(const int sampleNo, const int dataPointNo, const double value) {
  //
  // Get the number of samples and data-points per sample.
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int dataPointRank = getPointDataView().getRank();
  ShapeType dataPointShape = getPointDataView().getShape();
  if (numSamples*numDataPointsPerSample > 0) {
     //TODO: global error handling
     if ((sampleNo >= numSamples) || (sampleNo < 0 )) {
          throw DataException("Error - DataExpanded::copyDataPoint invalid sampleNo.");
     }
     if ((dataPointNo >= numDataPointsPerSample) || (dataPointNo < 0)) {
           throw DataException("Error - DataExpanded::copyDataPoint invalid dataPointNoInSample.");
     }
     DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNo);
     if (dataPointRank==0) {
         dataPointView()=value;
     } else if (dataPointRank==1) {
        for (int i=0; i<dataPointShape[0]; i++) {
            dataPointView(i)=value;
        }
     } else if (dataPointRank==2) {
        for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
              dataPointView(i,j)=value;
           }
        }
     } else if (dataPointRank==3) {
        for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
              for (int k=0; k<dataPointShape[2]; k++) {
                 dataPointView(i,j,k)=value;
              }
           }
        }
     } else if (dataPointRank==4) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
               for (int l=0; l<dataPointShape[3]; l++) {
                  dataPointView(i,j,k,l)=value;
               }
             }
           }
         }
     }
  }
}
void
DataExpanded::copyToDataPoint(const int sampleNo, const int dataPointNo, const boost::python::numeric::array& value) {
  //
  // Get the number of samples and data-points per sample.
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int dataPointRank = getPointDataView().getRank();
  ShapeType dataPointShape = getPointDataView().getShape();
  //
  // check rank:
  if (value.getrank()!=dataPointRank)
       throw DataException("Rank of numarray does not match Data object rank");
  if (numSamples*numDataPointsPerSample > 0) {
     //TODO: global error handling
     if ((sampleNo >= numSamples) || (sampleNo < 0 )) {
          throw DataException("Error - DataExpanded::copyDataPoint invalid sampleNo.");
     }
     if ((dataPointNo >= numDataPointsPerSample) || (dataPointNo < 0)) {
           throw DataException("Error - DataExpanded::copyDataPoint invalid dataPointNoInSample.");
     }
     DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNo);
     if (dataPointRank==0) {
         dataPointView()=extract<double>(value[0]);
     } else if (dataPointRank==1) {
        for (int i=0; i<dataPointShape[0]; i++) {
            dataPointView(i)=extract<double>(value[i]);
        }
     } else if (dataPointRank==2) {
        for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
              dataPointView(i,j)=extract<double>(value[i][j]);
           }
        }
     } else if (dataPointRank==3) {
        for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
              for (int k=0; k<dataPointShape[2]; k++) {
                 dataPointView(i,j,k)=extract<double>(value[i][j][k]);
              }
           }
        }
     } else if (dataPointRank==4) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
               for (int l=0; l<dataPointShape[3]; l++) {
                  dataPointView(i,j,k,l)=extract<double>(value[i][j][k][l]);
               }
             }
           }
         }
     }
  }
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
DataExpanded::symmetric(DataAbstract* ev)
{
  int sampleNo,dataPointNo;
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::symmetric: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::symmetric(thisView,getPointOffset(sampleNo,dataPointNo),
                                    evView,ev->getPointOffset(sampleNo,dataPointNo));
    }
  }
}
void
DataExpanded::nonsymmetric(DataAbstract* ev)
{
  int sampleNo,dataPointNo;
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::nonsymmetric: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::nonsymmetric(thisView,getPointOffset(sampleNo,dataPointNo),
                                    evView,ev->getPointOffset(sampleNo,dataPointNo));
    }
  }
}
void
DataExpanded::trace(DataAbstract* ev, int axis_offset)
{
  int sampleNo,dataPointNo;
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::trace: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::trace(thisView,getPointOffset(sampleNo,dataPointNo),
                                    evView,ev->getPointOffset(sampleNo,dataPointNo),axis_offset);
    }
  }
}

void
DataExpanded::transpose(DataAbstract* ev, int axis_offset)
{
  int sampleNo,dataPointNo;
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::transpose: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::transpose(thisView,getPointOffset(sampleNo,dataPointNo),
                                    evView,ev->getPointOffset(sampleNo,dataPointNo),axis_offset);
    }
  }
}

void
DataExpanded::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
  int sampleNo,dataPointNo;
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::swapaxes: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::swapaxes(thisView,getPointOffset(sampleNo,dataPointNo),
                                    evView,ev->getPointOffset(sampleNo,dataPointNo),axis0,axis1);
    }
  }
}
void
DataExpanded::eigenvalues(DataAbstract* ev)
{
  int sampleNo,dataPointNo;
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataExpanded::eigenvalues: casting to DataExpanded failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
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
  int sampleNo,dataPointNo;
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
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataArrayView::eigenvalues_and_eigenvectors(thisView,getPointOffset(sampleNo,dataPointNo),
                                                     evView,ev->getPointOffset(sampleNo,dataPointNo),
                                                     VView,V->getPointOffset(sampleNo,dataPointNo),
                                                     tol);
    }
  }
}

void
DataExpanded::setToZero(){
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataArrayView& thisView=getPointDataView();
  DataArrayView::ValueType::size_type n = thisView.noValues();
  double* p;
  int  sampleNo,dataPointNo, i;
  #pragma omp parallel for private(sampleNo,dataPointNo,p,i) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
        p=&(m_data[getPointOffset(sampleNo,dataPointNo)]);
        for (i=0; i<n ;++i) p[i]=0.;
    }
  }
}


void
DataExpanded::dump(const std::string fileName) const
{
   #ifdef PASO_MPI
   throw DataException("Error - DataExpanded:: dump is not implemented for MPI yet.");
   #endif
   #ifdef USE_NETCDF
   const int ldims=2+DataArrayView::maxRank;
   const NcDim* ncdims[ldims];
   NcVar *var, *ids;
   int rank = getPointDataView().getRank();
   int type=  getFunctionSpace().getTypeCode();
   int ndims =0;
   long dims[ldims];
   const double* d_ptr=&(m_data[0]);
   DataArrayView::ShapeType shape = getPointDataView().getShape();

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   NcFile dataFile(fileName.c_str(), NcFile::Replace);
   // check if writing was successful
   if (!dataFile.is_valid())
        throw DataException("Error - DataExpanded:: opening of netCDF file for output failed.");
   if (!dataFile.add_att("type_id",2) )
        throw DataException("Error - DataExpanded:: appending data type to netCDF file failed.");
   if (!dataFile.add_att("rank",rank) )
        throw DataException("Error - DataExpanded:: appending rank attribute to netCDF file failed.");
   if (!dataFile.add_att("function_space_type",type))
        throw DataException("Error - DataExpanded:: appending function space attribute to netCDF file failed.");
   ndims=rank+2;
   if ( rank >0 ) {
       dims[0]=shape[0];
       if (! (ncdims[0] = dataFile.add_dim("d0",shape[0])) )
            throw DataException("Error - DataExpanded:: appending ncdimsion 0 to netCDF file failed.");
   }
   if ( rank >1 ) {
       dims[1]=shape[1];
       if (! (ncdims[1] = dataFile.add_dim("d1",shape[1])) )
            throw DataException("Error - DataExpanded:: appending ncdimsion 1 to netCDF file failed.");
   }
   if ( rank >2 ) {
       dims[2]=shape[2];
       if (! (ncdims[2] = dataFile.add_dim("d2", shape[2])) )
            throw DataException("Error - DataExpanded:: appending ncdimsion 2 to netCDF file failed.");
   }
   if ( rank >3 ) {
       dims[3]=shape[3];
       if (! (ncdims[3] = dataFile.add_dim("d3", shape[3])) )
            throw DataException("Error - DataExpanded:: appending ncdimsion 3 to netCDF file failed.");
   }
   dims[rank]=getFunctionSpace().getNumDataPointsPerSample();
   if (! (ncdims[rank] = dataFile.add_dim("num_data_points_per_sample", dims[rank])) )
            throw DataException("Error - DataExpanded:: appending num_data_points_per_sample to netCDF file failed.");
   dims[rank+1]=getFunctionSpace().getNumSamples();
   if (! (ncdims[rank+1] = dataFile.add_dim("num_samples", dims[rank+1])) )
            throw DataException("Error - DataExpanded:: appending num_sample to netCDF file failed.");

   if (! ( ids = dataFile.add_var("id", ncInt, ncdims[rank+1])) )
        throw DataException("Error - DataExpanded:: appending reference id to netCDF file failed.");
   const int* ids_p=getFunctionSpace().borrowSampleReferenceIDs();
   if (! (ids->put(ids_p,dims[rank+1])) )
        throw DataException("Error - DataExpanded:: copy reference id  to netCDF buffer failed.");

   if (! ( var = dataFile.add_var("data", ncDouble, ndims, ncdims)) )
        throw DataException("Error - DataExpanded:: appending variable to netCDF file failed.");
   if (! (var->put(d_ptr,dims)) )
        throw DataException("Error - DataExpanded:: copy data to netCDF buffer failed.");
   #else
   throw DataException("Error - DataExpanded:: dump is not configured with netCDF. Please contact your installation manager.");
   #endif
}

void
DataExpanded::setTaggedValue(int tagKey,
                             const DataArrayView& value)
{
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int sampleNo,dataPointNo, i;
  DataArrayView& thisView=getPointDataView();
  DataArrayView::ValueType::size_type n = thisView.noValues();
  double* p,*in=&(value.getData()[0]);
  
  if (value.noValues() != n) {
    throw DataException("Error - DataExpanded::setTaggedValue: number of input values does not match number of values per data points.");
  }

  #pragma omp parallel for private(sampleNo,dataPointNo,p,i) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    if (getFunctionSpace().getTagFromSampleNo(sampleNo) == tagKey ) {
        for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            p=&(m_data[getPointOffset(sampleNo,dataPointNo)]);
            for (i=0; i<n ;++i) p[i]=in[i];
        }
    }
  }
}


}  // end of namespace