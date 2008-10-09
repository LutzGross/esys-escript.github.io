
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "Data.h"
#include "DataExpanded.h"
#include "DataException.h"
#include "DataConstant.h"
#include "DataTagged.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif

#include <boost/python/extract.hpp>
#include "DataMaths.h"

using namespace std;
using namespace boost::python;
using namespace boost;
using namespace escript::DataTypes;

namespace escript {

DataExpanded::DataExpanded(const boost::python::numeric::array& value,
                           const FunctionSpace& what)
  : parent(what,DataTypes::shapeFromNumArray(value))
{
  //
  // initialise the data array for this object
  initialise(what.getNumSamples(),what.getNumDPPSample());
  //
  // copy the given value to every data point
  copy(value);
}

DataExpanded::DataExpanded(const DataExpanded& other)
  : parent(other.getFunctionSpace(), other.getShape()),
  m_data(other.m_data)
{
}

DataExpanded::DataExpanded(const DataConstant& other)
  : parent(other.getFunctionSpace(), other.getShape())
{
  //
  // initialise the data array for this object
  initialise(other.getNumSamples(),other.getNumDPPSample());
  //
  // DataConstant only has one value, copy this to every data point
  copy(other);
}

DataExpanded::DataExpanded(const DataTagged& other)
  : parent(other.getFunctionSpace(), other.getShape())
{
  //
  // initialise the data array for this object
  initialise(other.getNumSamples(),other.getNumDPPSample());
  //
  // for each data point in this object, extract and copy the corresponding data
  // value from the given DataTagged object
  int i,j;
  DataTypes::ValueType::size_type numRows=m_data.getNumRows();
  DataTypes::ValueType::size_type numCols=m_data.getNumCols();
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numRows;i++) {
    for (j=0;j<numCols;j++) {
      try {
           DataTypes::copyPoint(getVector(), getPointOffset(i,j), getNoValues(),
                                other.getVector(),
                                other.getPointOffset(i,j));
      }
      catch (std::exception& e) {
        cout << e.what() << endl;
      }
    }
  }
}

DataExpanded::DataExpanded(const DataExpanded& other,
                           const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
  //
  // get the shape of the slice
//   DataTypes::ShapeType shape(DataTypes::getResultSliceShape(region));
  //
  // initialise this Data object to the shape of the slice
  initialise(other.getNumSamples(),other.getNumDPPSample());
  //
  // copy the data
  DataTypes::RegionLoopRangeType region_loop_range=DataTypes::getSliceRegionLoopRange(region);
  DataTypes::ValueType::size_type numRows=m_data.getNumRows();
  DataTypes::ValueType::size_type numCols=m_data.getNumCols();
  int i,j;
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numRows;i++) {
    for (j=0;j<numCols;j++) {
      try {
//         getPointDataView().copySlice(getPointOffset(i,j),
//                                      other.getPointDataView(),
//                                      other.getPointOffset(i,j),
//                                      region_loop_range);
        DataTypes::copySlice(getVector(),getShape(),getPointOffset(i,j),
                                     other.getVector(),
				     other.getShape(),
                                     other.getPointOffset(i,j),
                                     region_loop_range);
      }
      catch (std::exception& e) {
        cout << e.what() << endl;
      }
    }
  }
}

// DataExpanded::DataExpanded(const DataArrayView& value,
//                            const FunctionSpace& what)
//   : DataAbstract(what)
// {
//   //
//   // get the shape of the given data value
//   DataTypes::ShapeType tempShape=value.getShape();
//   //
//   // initialise this Data object to the shape of the given data value
//   initialise(tempShape,what.getNumSamples(),what.getNumDPPSample());
//   //
//   // copy the given value to every data point
//   copy(value);
// }

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::ValueType &data)
  : parent(what,shape)
{
  EsysAssert(data.size()%getNoValues()==0,
                 "DataExpanded Constructor - size of supplied data is not a multiple of shape size.");

  if (data.size()==getNoValues())
  {
     ValueType& vec=m_data.getData();
     //
     // create the view of the data
     initialise(what.getNumSamples(),what.getNumDPPSample());
     // now we copy this value to all elements
     for (int i=0;i<getLength();)
     {
	for (int j=0;j<getNoValues();++j,++i)
	{
	    vec[i]=data[j];
	}
     }
  }
  else
  {
     //
     // copy the data in the correct format
     m_data.getData()=data;
  }


}

DataExpanded::~DataExpanded()
{
}

DataAbstract*
DataExpanded::deepCopy()
{
  return new DataExpanded(*this);
}


DataAbstract*
DataExpanded::getSlice(const DataTypes::RegionType& region) const
{
  return new DataExpanded(*this,region);
}

void
DataExpanded::setSlice(const DataAbstract* value,
                       const DataTypes::RegionType& region)
{
  const DataExpanded* tempDataExp=dynamic_cast<const DataExpanded*>(value);
  if (tempDataExp==0) {
    throw DataException("Programming error - casting to DataExpanded.");
  }
  //
  // get shape of slice
  DataTypes::ShapeType shape(DataTypes::getResultSliceShape(region));
  DataTypes::RegionLoopRangeType region_loop_range=DataTypes::getSliceRegionLoopRange(region);
  //
  // check shape
  if (getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (tempDataExp->getRank()>0 && !DataTypes::checkShape(value->getShape(), shape)) {
    throw DataException (DataTypes::createShapeErrorMessage(
		"Error - Couldn't copy slice due to shape mismatch.",shape, value->getShape()));
  }
  //
  // copy the data from the slice into this object
  DataTypes::ValueType::size_type numRows=m_data.getNumRows();
  DataTypes::ValueType::size_type numCols=m_data.getNumCols();
  int i, j;
  ValueType& vec=getVector();
  const ShapeType& mshape=getShape();
  const ValueType& tVec=tempDataExp->getVector();
  const ShapeType& tShape=tempDataExp->getShape();
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numRows;i++) {
    for (j=0;j<numCols;j++) {
/*      getPointDataView().copySliceFrom(getPointOffset(i,j),
                                       tempDataExp->getPointDataView(),
                                       tempDataExp->getPointOffset(i,j),
                                       region_loop_range);*/
        DataTypes::copySliceFrom(vec,mshape,getPointOffset(i,j),
                                       tVec,
				       tShape,
                                       tempDataExp->getPointOffset(i,j),
                                       region_loop_range);

    }
  }
}

void
DataExpanded::copy(const DataConstant& value)
{
  EsysAssert((checkShape(getShape(), value.getShape())),
                 createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",value.getShape(),getShape()));

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
      //getPointDataView().copy(getPointOffset(i,j),value);
      DataTypes::copyPoint(getVector(), getPointOffset(i,j), getNoValues(), value.getVector(), 0);
    }
  }
}


// void
// DataExpanded::copy(const DataArrayView& value)
// {
//   //
//   // copy a single value to every data point in this object
//   int nRows=m_data.getNumRows();
//   int nCols=m_data.getNumCols();
//   int i,j;
//   #pragma omp parallel for private(i,j) schedule(static)
//   for (i=0;i<nRows;i++) {
//     for (j=0;j<nCols;j++) {
//       // NOTE: An exception may be thown from this call if
//       // DOASSERT is on which of course will play
//       // havoc with the omp threads. Run single threaded
//       // if using DOASSERT.
//       getPointDataView().copy(getPointOffset(i,j),value);
//     }
//   }
// }

void
DataExpanded::copy(const boost::python::numeric::array& value)
{

  // extract the shape of the numarray
  DataTypes::ShapeType tempShape;
  for (int i=0; i < value.getrank(); i++) {
    tempShape.push_back(extract<int>(value.getshape()[i]));
  }

  // get the space for the data vector
//   int len = DataTypes::noValues(tempShape);
//   DataVector temp_data(len, 0.0, len);
//   DataArrayView temp_dataView(temp_data, tempShape);
//   temp_dataView.copy(value);

  //
  // check the input shape matches this shape
  if (!DataTypes::checkShape(getShape(),tempShape)) {
    throw DataException(DataTypes::createShapeErrorMessage(
                        "Error - (DataExpanded) Cannot copy due to shape mismatch.",
                        tempShape,getShape()));
  }
  //
  // now copy over the data
  //copy(temp_dataView);
  getVector().copyFromNumArray(value);
}


void
DataExpanded::initialise(int noSamples,
                         int noDataPointsPerSample)
{
  if (noSamples==0)		//retain the default empty object
  {
     return;
  }
  //
  // resize data array to the required size
  m_data.resize(noSamples,noDataPointsPerSample,getNoValues());
}

string
DataExpanded::toString() const
{
  stringstream temp;
  FunctionSpace fs=getFunctionSpace();

  int offset=0;
  for (int i=0;i<m_data.getNumRows();i++) {
    for (int j=0;j<m_data.getNumCols();j++) {
      offset=getPointOffset(i,j);
      stringstream suffix;
      suffix << "( id: " << i << ", ref: " << fs.getReferenceIDOfSample(i) << ", pnt: " << j << ")";
      temp << DataTypes::pointToString(getVector(),getShape(),offset,suffix.str());
      if (!(i==(m_data.getNumRows()-1) && j==(m_data.getNumCols()-1))) {
        temp << endl;
      }
    }
  }
  string result=temp.str();
  if (result.empty())
  {
    return "(data contains no samples)\n";
  }
  return temp.str();
}

DataTypes::ValueType::size_type
DataExpanded::getPointOffset(int sampleNo,
                             int dataPointNo) const
{
  return m_data.index(sampleNo,dataPointNo);
}

DataTypes::ValueType::size_type
DataExpanded::getLength() const
{
  return m_data.size();
}



void
DataExpanded::copyToDataPoint(const int sampleNo, const int dataPointNo, const double value) {
  //
  // Get the number of samples and data-points per sample.
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int dataPointRank = getRank();
  ShapeType dataPointShape = getShape();
  if (numSamples*numDataPointsPerSample > 0) {
     //TODO: global error handling
     if ((sampleNo >= numSamples) || (sampleNo < 0 )) {
          throw DataException("Error - DataExpanded::copyDataPoint invalid sampleNo.");
     }
     if ((dataPointNo >= numDataPointsPerSample) || (dataPointNo < 0)) {
           throw DataException("Error - DataExpanded::copyDataPoint invalid dataPointNoInSample.");
     }
     ValueType::size_type offset = getPointOffset(sampleNo, dataPointNo);
     ValueType& vec=getVector();
     if (dataPointRank==0) {
         vec[offset]=value;
     } else if (dataPointRank==1) {
        for (int i=0; i<dataPointShape[0]; i++) {
            vec[offset+i]=value;
        }
     } else if (dataPointRank==2) {
        for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
              vec[offset+getRelIndex(dataPointShape,i,j)]=value;
           }
        }
     } else if (dataPointRank==3) {
        for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
              for (int k=0; k<dataPointShape[2]; k++) {
                 vec[offset+getRelIndex(dataPointShape,i,j,k)]=value;
              }
           }
        }
     } else if (dataPointRank==4) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
               for (int l=0; l<dataPointShape[3]; l++) {
                  vec[offset+getRelIndex(dataPointShape,i,j,k,l)]=value;
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
  int dataPointRank = getRank();
  const ShapeType& shape = getShape();
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
     ValueType& vec=getVector();
     if (dataPointRank==0) {
         vec[0]=extract<double>(value[0]);
     } else if (dataPointRank==1) {
        for (int i=0; i<shape[0]; i++) {
            vec[i]=extract<double>(value[i]);
        }
     } else if (dataPointRank==2) {
        for (int i=0; i<shape[0]; i++) {
           for (int j=0; j<shape[1]; j++) {
              vec[getRelIndex(shape,i,j)]=extract<double>(value[i][j]);
           }
        }
     } else if (dataPointRank==3) {
        for (int i=0; i<shape[0]; i++) {
           for (int j=0; j<shape[1]; j++) {
              for (int k=0; k<shape[2]; k++) {
                 vec[getRelIndex(shape,i,j,k)]=extract<double>(value[i][j][k]);
              }
           }
        }
     } else if (dataPointRank==4) {
         for (int i=0; i<shape[0]; i++) {
           for (int j=0; j<shape[1]; j++) {
             for (int k=0; k<shape[2]; k++) {
               for (int l=0; l<shape[3]; l++) {
                  vec[getRelIndex(shape,i,j,k,l)]=extract<double>(value[i][j][k][l]);
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
  int dataPointRank = getRank();
  const ShapeType& dataPointShape = getShape();
  //
  // check rank:
  if (value.getrank()!=dataPointRank+1)
       throw DataException("Rank of numarray does not match Data object rank");
  if (value.getshape()[0]!=numSamples*numDataPointsPerSample)
       throw DataException("leading dimension of numarray is too small");
  //
  ValueType& vec=getVector();
  int dataPoint = 0;
  for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
      ValueType::size_type offset=getPointOffset(sampleNo, dataPointNo);
      if (dataPointRank==0) {
         vec[offset]=extract<double>(value[dataPoint]);
      } else if (dataPointRank==1) {
         for (int i=0; i<dataPointShape[0]; i++) {
            vec[offset+i]=extract<double>(value[dataPoint][i]);
         }
      } else if (dataPointRank==2) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
	     vec[offset+getRelIndex(dataPointShape,i,j)]=extract<double>(value[dataPoint][i][j]);
           }
         }
       } else if (dataPointRank==3) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
		 vec[offset+getRelIndex(dataPointShape,i,j,k)]=extract<double>(value[dataPoint][i][j][k]);
             }
           }
         }
       } else if (dataPointRank==4) {
         for (int i=0; i<dataPointShape[0]; i++) {
           for (int j=0; j<dataPointShape[1]; j++) {
             for (int k=0; k<dataPointShape[2]; k++) {
               for (int l=0; l<dataPointShape[3]; l++) {
                 vec[offset+getRelIndex(dataPointShape,i,j,k,l)]=extract<double>(value[dataPoint][i][j][k][l]);
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::symmetric(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo));
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::nonsymmetric(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo));
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::trace(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo),axis_offset);
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::transpose(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo),axis_offset);
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::swapaxes(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo),axis0,axis1);
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::eigenvalues(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo));
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
  ValueType& vec=getVector();
  const ShapeType& shape=getShape();
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  ValueType& VVec=temp_V->getVector();
  const ShapeType& VShape=temp_V->getShape();
  #pragma omp parallel for private(sampleNo,dataPointNo) schedule(static)
  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
    for (dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
         DataMaths::eigenvalues_and_eigenvectors(vec,shape,getPointOffset(sampleNo,dataPointNo),
                                    evVec,evShape,ev->getPointOffset(sampleNo,dataPointNo),
                                    VVec, VShape,V->getPointOffset(sampleNo,dataPointNo),tol);
    }
  }
}

void
DataExpanded::setToZero(){
// TODO: Surely there is a more efficient way to do this????
// Why is there no memset here? Parallel issues?
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  DataTypes::ValueType::size_type n = getNoValues();
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

/* Append MPI rank to file name if multiple MPI processes */
char *Escript_MPI_appendRankToFileName(const char *fileName, int mpi_size, int mpi_rank) {
  /* Make plenty of room for the mpi_rank number and terminating '\0' */
  char *newFileName = (char *)malloc(strlen(fileName)+20);
  strncpy(newFileName, fileName, strlen(fileName)+1);
  if (mpi_size>1) sprintf(newFileName+strlen(newFileName), ".%04d", mpi_rank);
  return(newFileName);
}

void
DataExpanded::dump(const std::string fileName) const
{
   #ifdef USE_NETCDF
   const int ldims=2+DataTypes::maxRank;
   const NcDim* ncdims[ldims];
   NcVar *var, *ids;
   int rank = getRank();
   int type=  getFunctionSpace().getTypeCode();
   int ndims =0;
   long dims[ldims];
   const double* d_ptr=&(m_data[0]);
   const DataTypes::ShapeType& shape = getShape();
   int mpi_iam=getFunctionSpace().getDomain()->getMPIRank();
   int mpi_num=getFunctionSpace().getDomain()->getMPISize();
#ifdef PASO_MPI
   MPI_Status status;
#endif

#ifdef PASO_MPI
   /* Serialize NetCDF I/O */
   if (mpi_iam>0) MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81801, MPI_COMM_WORLD, &status);
#endif

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   char *newFileName = Escript_MPI_appendRankToFileName(fileName.c_str(), mpi_num, mpi_iam);
   NcFile dataFile(newFileName, NcFile::Replace);
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
            throw DataException("Error - DataExpanded:: appending ncdimension 0 to netCDF file failed.");
   }
   if ( rank >1 ) {
       dims[1]=shape[1];
       if (! (ncdims[1] = dataFile.add_dim("d1",shape[1])) )
            throw DataException("Error - DataExpanded:: appending ncdimension 1 to netCDF file failed.");
   }
   if ( rank >2 ) {
       dims[2]=shape[2];
       if (! (ncdims[2] = dataFile.add_dim("d2", shape[2])) )
            throw DataException("Error - DataExpanded:: appending ncdimension 2 to netCDF file failed.");
   }
   if ( rank >3 ) {
       dims[3]=shape[3];
       if (! (ncdims[3] = dataFile.add_dim("d3", shape[3])) )
            throw DataException("Error - DataExpanded:: appending ncdimension 3 to netCDF file failed.");
   }
   dims[rank]=getFunctionSpace().getNumDataPointsPerSample();
   if (! (ncdims[rank] = dataFile.add_dim("num_data_points_per_sample", dims[rank])) )
            throw DataException("Error - DataExpanded:: appending num_data_points_per_sample to netCDF file failed.");
   dims[rank+1]=getFunctionSpace().getNumSamples();
   if (! (ncdims[rank+1] = dataFile.add_dim("num_samples", dims[rank+1])) )
            throw DataException("Error - DataExpanded:: appending num_sample to netCDF file failed.");

   if (getFunctionSpace().getNumSamples()>0)
   {

     if (! ( ids = dataFile.add_var("id", ncInt, ncdims[rank+1])) )
        throw DataException("Error - DataExpanded:: appending reference id to netCDF file failed.");
     const int* ids_p=getFunctionSpace().borrowSampleReferenceIDs();
     if (! (ids->put(ids_p,dims[rank+1])) )
        throw DataException("Error - DataExpanded:: copy reference id  to netCDF buffer failed.");
     if (! ( var = dataFile.add_var("data", ncDouble, ndims, ncdims)) )
        throw DataException("Error - DataExpanded:: appending variable to netCDF file failed.");
     if (! (var->put(d_ptr,dims)) )
        throw DataException("Error - DataExpanded:: copy data to netCDF buffer failed.");
   }
#ifdef PASO_MPI
   if (mpi_iam<mpi_num-1) MPI_Send(&ndims, 0, MPI_INT, mpi_iam+1, 81801, MPI_COMM_WORLD);
#endif
   #else
   throw DataException("Error - DataExpanded:: dump is not configured with netCDF. Please contact your installation manager.");
   #endif
}

void  
DataExpanded::setTaggedValue(int tagKey,
 	       const DataTypes::ShapeType& pointshape,
               const DataTypes::ValueType& value,
	       int dataOffset)
{
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int sampleNo,dataPointNo, i;
  DataTypes::ValueType::size_type n = getNoValues();
  double* p;
  const double* in=&value[0+dataOffset];
  
  if (value.size() != n) {
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


void
DataExpanded::reorderByReferenceIDs(int *reference_ids)
{
  int numSamples = getNumSamples();
  DataTypes::ValueType::size_type n = getNoValues() * getNumDPPSample();
  int sampleNo, sampleNo2,i;
  double* p,*p2;
  register double rtmp;
  FunctionSpace fs=getFunctionSpace();

  for (sampleNo = 0; sampleNo < numSamples; sampleNo++) {
     const int id_in=reference_ids[sampleNo];
     const int id=fs.getReferenceIDOfSample(sampleNo);
     if (id!=id_in) {
         bool matched=false;
         for (sampleNo2 = sampleNo+1; sampleNo2 < numSamples; sampleNo2++) {
              if (id == reference_ids[sampleNo2]) {
                 p=&(m_data[getPointOffset(sampleNo,0)]);
                 p2=&(m_data[getPointOffset(sampleNo2,0)]);
                 for (i=0; i<n ;i++) {
                         rtmp=p[i];
                         p[i]=p2[i];
                         p2[i]=rtmp;
                 }
                 reference_ids[sampleNo]=id;
                 reference_ids[sampleNo2]=id_in;
                 matched=true;
                 break;
              }
         }
         if (! matched) {
            throw DataException("Error - DataExpanded::reorderByReferenceIDs: unable to reorder sample data by reference ids");
         }
     }
   }
}

DataTypes::ValueType&
DataExpanded::getVector()
{
	return m_data.getData();
}

const DataTypes::ValueType&
DataExpanded::getVector() const
{
	return m_data.getData();
}

}  // end of namespace
