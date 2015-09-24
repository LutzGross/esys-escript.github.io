
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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

#define ESNEEDPYTHON
#include <esysUtils/first.h>
#include <esysUtils/Esys_MPI.h>

#include "Data.h"
#include "DataConstant.h"
#include "DataException.h"
#include "DataExpanded.h"
#include "DataMaths.h"
#include "DataTagged.h"

#include <limits>
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

using namespace std;
using namespace escript::DataTypes;

#define CHECK_FOR_EX_WRITE do {\
    if (!checkNoSharing()) {\
        std::ostringstream ss;\
        ss << "Attempt to modify shared object. Line " << __LINE__ << " in "\
           << __FILE__;\
        abort();\
        throw DataException(ss.str());\
    }\
} while(0)


namespace escript {

DataExpanded::DataExpanded(const WrappedArray& value,
                           const FunctionSpace& what)
  : parent(what,value.getShape())
{
    // initialise the data array for this object
    initialise(what.getNumSamples(),what.getNumDPPSample());
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
    // initialise the data array for this object
    initialise(other.getNumSamples(),other.getNumDPPSample());
    // DataConstant only has one value, copy this to every data point
    copy(other);
}

DataExpanded::DataExpanded(const DataTagged& other)
  : parent(other.getFunctionSpace(), other.getShape())
{
    // initialise the data array for this object
    initialise(other.getNumSamples(),other.getNumDPPSample());
    // for each data point in this object, extract and copy the corresponding
    // data value from the given DataTagged object
    DataTypes::ValueType::size_type numRows=m_data.getNumRows();
    DataTypes::ValueType::size_type numCols=m_data.getNumCols();
#pragma omp parallel for
    for (int i=0; i<numRows; i++) {
        for (int j=0; j<numCols; j++) {
            try {
                DataTypes::copyPoint(getVectorRW(), getPointOffset(i,j),
                                     getNoValues(), other.getVectorRO(),
                                     other.getPointOffset(i,j));
            } catch (std::exception& e) {
                cerr << e.what() << endl;
            }
        }
    }
}

DataExpanded::DataExpanded(const DataExpanded& other,
                           const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
    // initialise this Data object to the shape of the slice
    initialise(other.getNumSamples(),other.getNumDPPSample());
    // copy the data
    DataTypes::RegionLoopRangeType region_loop_range =
                                    DataTypes::getSliceRegionLoopRange(region);
    DataTypes::ValueType::size_type numRows=m_data.getNumRows();
    DataTypes::ValueType::size_type numCols=m_data.getNumCols();
#pragma omp parallel for
    for (int i=0; i<numRows; i++) {
        for (int j=0; j<numCols; j++) {
            try {
                DataTypes::copySlice(getVectorRW(), getShape(),
                                     getPointOffset(i,j), other.getVectorRO(),
                                     other.getShape(),
                                     other.getPointOffset(i,j),
                                     region_loop_range);
            } catch (std::exception& e) {
                cerr << e.what() << endl;
            }
        }
    }
}

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::ValueType &data)
  : parent(what,shape)
{
    EsysAssert(data.size()%getNoValues()==0,
                 "DataExpanded Constructor - size of supplied data is not a multiple of shape size.");

    if (data.size() == getNoValues()) {
        ValueType& vec=m_data.getData();
        // create the view of the data
        initialise(what.getNumSamples(),what.getNumDPPSample());
        // now we copy this value to all elements
        for (int i=0; i<getLength();) {
            for (unsigned int j=0;j<getNoValues();++j,++i) {
                vec[i]=data[j];
            }
        }
    } else {
        // copy the data in the correct format
        m_data.getData() = data;
    }
}

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const double v)
  : parent(what,shape)
{
    ValueType& vec=m_data.getData();
    // create the view of the data
    initialise(what.getNumSamples(),what.getNumDPPSample());
    // now we copy this value to all elements
    const int L=getLength();
#pragma omp parallel for
    for (int i=0; i<L; ++i) {
        vec[i]=v;
    }
}

DataExpanded::~DataExpanded()
{
}

DataAbstract* DataExpanded::deepCopy()
{
    return new DataExpanded(*this);
}

DataAbstract* DataExpanded::getSlice(const DataTypes::RegionType& region) const
{
    return new DataExpanded(*this,region);
}

void DataExpanded::setSlice(const DataAbstract* value,
                            const DataTypes::RegionType& region)
{
    const DataExpanded* tempDataExp=dynamic_cast<const DataExpanded*>(value);
    if (!tempDataExp)
        throw DataException("Programming error - casting to DataExpanded.");

    CHECK_FOR_EX_WRITE;

    // get shape of slice
    DataTypes::ShapeType shape(DataTypes::getResultSliceShape(region));
    DataTypes::RegionLoopRangeType region_loop_range =
                                 DataTypes::getSliceRegionLoopRange(region);
    // check shape
    if (getRank() != region.size())
        throw DataException("Error - Invalid slice region.");

    if (tempDataExp->getRank() > 0 && !DataTypes::checkShape(value->getShape(), shape))
        throw DataException(DataTypes::createShapeErrorMessage(
                "Error - Couldn't copy slice due to shape mismatch.",
                shape, value->getShape()));

    // copy the data from the slice into this object
    DataTypes::ValueType::size_type numRows = m_data.getNumRows();
    DataTypes::ValueType::size_type numCols = m_data.getNumCols();
    ValueType& vec=getVectorRW();
    const ShapeType& mshape=getShape();
    const ValueType& tVec=tempDataExp->getVectorRO();
    const ShapeType& tShape=tempDataExp->getShape();
#pragma omp parallel for
    for (int i=0; i<numRows; i++) {
        for (int j=0; j<numCols; j++) {
            DataTypes::copySliceFrom(vec, mshape, getPointOffset(i,j), tVec,
                                     tShape, tempDataExp->getPointOffset(i,j),
                                     region_loop_range);
        }
    }
}

void DataExpanded::copy(const DataConstant& value)
{
    EsysAssert((checkShape(getShape(), value.getShape())),
                 createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.", value.getShape(), getShape()));

    // copy a single value to every data point in this object
    int nRows=m_data.getNumRows();
    int nCols=m_data.getNumCols();
#pragma omp parallel for
    for (int i=0; i<nRows; i++) {
        for (int j=0; j<nCols; j++) {
            DataTypes::copyPoint(getVectorRW(), getPointOffset(i,j),
                                 getNoValues(), value.getVectorRO(), 0);
        }
    }
}

void DataExpanded::copy(const WrappedArray& value)
{
    // check the input shape matches this shape
    if (!DataTypes::checkShape(getShape(),value.getShape()))
        throw DataException(DataTypes::createShapeErrorMessage(
                   "Error - (DataExpanded) Cannot copy due to shape mismatch.",
                   value.getShape(),getShape()));
    getVectorRW().copyFromArray(value, getNumDPPSample()*getNumSamples());
}

void DataExpanded::initialise(int noSamples, int noDataPointsPerSample)
{
    if (noSamples==0) //retain the default empty object
        return;

    // resize data array to the required size
    m_data.resize(noSamples, noDataPointsPerSample, getNoValues());
}

bool DataExpanded::hasNaN() const
{
    bool haveNaN = false;
    const ValueType& v = m_data.getData();
#pragma omp parallel for
    for (ValueType::size_type i=0; i<v.size(); ++i) {
        if (nancheck(v[i])) {
#pragma omp critical
            {
                haveNaN=true;
            }
        }
    }
    return haveNaN;
}

void DataExpanded::replaceNaN(double value)
{
#pragma omp parallel for
    for (ValueType::size_type i=0; i<m_data.size(); ++i) {
        if (nancheck(m_data[i])) {
            m_data[i] = value;
        }
    }
}

string DataExpanded::toString() const
{
    stringstream ss;
    FunctionSpace fs=getFunctionSpace();

    int offset=0;
    for (int i=0; i<m_data.getNumRows(); i++) {
        for (int j=0; j<m_data.getNumCols(); j++) {
            offset = getPointOffset(i,j);
            stringstream suffix;
            suffix << "( id: " << i << ", ref: "
                   << fs.getReferenceIDOfSample(i) << ", pnt: " << j << ")";
            ss << DataTypes::pointToString(getVectorRO(), getShape(), offset,
                                           suffix.str());
            if (!(i==(m_data.getNumRows()-1) && j==(m_data.getNumCols()-1))) {
                ss << endl;
            }
        }
    }
    string result(ss.str());
    if (result.empty())
        return "(data contains no samples)\n";

    return result;
}

DataTypes::ValueType::size_type DataExpanded::getPointOffset(int sampleNo,
                                                        int dataPointNo) const
{
    return m_data.index(sampleNo,dataPointNo);
}

DataTypes::ValueType::size_type DataExpanded::getPointOffset(int sampleNo,
                                                             int dataPointNo)
{
    return m_data.index(sampleNo,dataPointNo);
}

DataTypes::ValueType::size_type DataExpanded::getLength() const
{
    return m_data.size();
}

void DataExpanded::copyToDataPoint(int sampleNo, int dataPointNo, double value)
{
    CHECK_FOR_EX_WRITE;
    // Get the number of samples and data-points per sample.
    int numSamples = getNumSamples();
    int numDataPointsPerSample = getNumDPPSample();
    int dataPointRank = getRank();
    ShapeType dataPointShape = getShape();
    if (numSamples*numDataPointsPerSample > 0) {
        //TODO: global error handling
        if (sampleNo >= numSamples || sampleNo < 0)
            throw DataException("DataExpanded::copyDataPoint: invalid sampleNo.");
        if (dataPointNo >= numDataPointsPerSample || dataPointNo < 0)
            throw DataException("DataExpanded::copyDataPoint: invalid dataPointNo.");

        ValueType::size_type offset = getPointOffset(sampleNo, dataPointNo);
        ValueType& vec = getVectorRW();
        if (dataPointRank==0) {
            vec[offset] = value;
        } else if (dataPointRank==1) {
            for (int i=0; i<dataPointShape[0]; i++) {
                vec[offset+i] = value;
            }
        } else if (dataPointRank==2) {
            for (int i=0; i<dataPointShape[0]; i++) {
                for (int j=0; j<dataPointShape[1]; j++) {
                    vec[offset+getRelIndex(dataPointShape,i,j)] = value;
                }
            }
        } else if (dataPointRank==3) {
            for (int i=0; i<dataPointShape[0]; i++) {
                for (int j=0; j<dataPointShape[1]; j++) {
                    for (int k=0; k<dataPointShape[2]; k++) {
                        vec[offset+getRelIndex(dataPointShape,i,j,k)] = value;
                    }
                }
            }
        } else if (dataPointRank==4) {
            for (int i=0; i<dataPointShape[0]; i++) {
                for (int j=0; j<dataPointShape[1]; j++) {
                    for (int k=0; k<dataPointShape[2]; k++) {
                        for (int l=0; l<dataPointShape[3]; l++) {
                            vec[offset+getRelIndex(dataPointShape,i,j,k,l)] = value;
                        }
                    }
                }
            }
        }
    }
}

void DataExpanded::copyToDataPoint(int sampleNo, int dataPointNo,
                                   const WrappedArray& value)
{
    CHECK_FOR_EX_WRITE;
    // Get the number of samples and data-points per sample
    int numSamples = getNumSamples();
    int numDataPointsPerSample = getNumDPPSample();
    // check rank
    if (value.getRank() != getRank())
        throw DataException("Rank of value does not match Data object rank");

    if (numSamples*numDataPointsPerSample > 0) {
        //TODO: global error handling
        if (sampleNo >= numSamples || sampleNo < 0)
            throw DataException("DataExpanded::copyDataPoint: invalid sampleNo.");
        if (dataPointNo >= numDataPointsPerSample || dataPointNo < 0)
            throw DataException("DataExpanded::copyDataPoint: invalid dataPointNoInSample.");

        ValueType::size_type offset = getPointOffset(sampleNo, dataPointNo);
        ValueType& vec = getVectorRW();
        vec.copyFromArrayToOffset(value, offset, 1);
    }
}

void DataExpanded::symmetric(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev = dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::symmetric: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec = getVectorRO();
    const ShapeType& shape = getShape();
    ValueType& evVec = temp_ev->getVectorRW();
    const ShapeType& evShape = temp_ev->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::symmetric(vec, shape,
                    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
                    ev->getPointOffset(sampleNo,dataPointNo));
        }
    }
}

void DataExpanded::nonsymmetric(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::nonsymmetric: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec = getVectorRO();
    const ShapeType& shape = getShape();
    ValueType& evVec = temp_ev->getVectorRW();
    const ShapeType& evShape = temp_ev->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::nonsymmetric(vec, shape,
                    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
                    ev->getPointOffset(sampleNo,dataPointNo));
        }
    }
}

void DataExpanded::trace(DataAbstract* ev, int axis_offset)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev = dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::trace: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec=getVectorRO();
    const ShapeType& shape=getShape();
    ValueType& evVec=temp_ev->getVectorRW();
    const ShapeType& evShape=temp_ev->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::trace(vec, shape, getPointOffset(sampleNo,dataPointNo),
                     evVec, evShape, ev->getPointOffset(sampleNo,dataPointNo),
                     axis_offset);
        }
    }
}

void DataExpanded::transpose(DataAbstract* ev, int axis_offset)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::transpose: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec = getVectorRO();
    const ShapeType& shape = getShape();
    ValueType& evVec = temp_ev->getVectorRW();
    const ShapeType& evShape = temp_ev->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::transpose(vec, shape,
                    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
                    ev->getPointOffset(sampleNo,dataPointNo), axis_offset);
        }
    }
}

void DataExpanded::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("Error - DataExpanded::swapaxes: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec=getVectorRO();
    const ShapeType& shape=getShape();
    ValueType& evVec=temp_ev->getVectorRW();
    const ShapeType& evShape=temp_ev->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::swapaxes(vec, shape,
                    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
                    ev->getPointOffset(sampleNo,dataPointNo), axis0, axis1);
        }
    }
}

void DataExpanded::eigenvalues(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::eigenvalues: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec=getVectorRO();
    const ShapeType& shape=getShape();
    ValueType& evVec=temp_ev->getVectorRW();
    const ShapeType& evShape=temp_ev->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::eigenvalues(vec, shape,
                    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
                    ev->getPointOffset(sampleNo,dataPointNo));
        }
    }
}

void DataExpanded::eigenvalues_and_eigenvectors(DataAbstract* ev,
                                                DataAbstract* V, double tol)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::eigenvalues_and_eigenvectors: casting to DataExpanded failed (probably a programming error).");

    DataExpanded* temp_V=dynamic_cast<DataExpanded*>(V);
    if (!temp_V)
        throw DataException("DataExpanded::eigenvalues_and_eigenvectors: casting to DataExpanded failed (probably a programming error).");

    const ValueType& vec = getVectorRO();
    const ShapeType& shape = getShape();
    ValueType& evVec = temp_ev->getVectorRW();
    const ShapeType& evShape = temp_ev->getShape();
    ValueType& VVec = temp_V->getVectorRW();
    const ShapeType& VShape = temp_V->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            DataMaths::eigenvalues_and_eigenvectors(vec, shape,
                    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
                    ev->getPointOffset(sampleNo,dataPointNo), VVec, VShape,
                    V->getPointOffset(sampleNo,dataPointNo), tol);
        }
    }
}


int DataExpanded::matrixInverse(DataAbstract* out) const
{
    DataExpanded* temp=dynamic_cast<DataExpanded*>(out);
    if (!temp)
        throw DataException("DataExpanded::matrixInverse: casting to DataExpanded failed (probably a programming error).");

    if (getRank() != 2)
        throw DataException("DataExpanded::matrixInverse: input must be rank 2.");

    const int numdpps=getNumDPPSample();
    const int numSamples = getNumSamples();
    const ValueType& vec=m_data.getData();
    int errcode=0;
#pragma omp parallel
    {
        int errorcode=0;
        LapackInverseHelper h(getShape()[0]);
#pragma omp for
        for (int sampleNo = 0; sampleNo < numSamples; sampleNo++)
        {
            // not sure I like all those virtual calls to getPointOffset
            DataTypes::ValueType::size_type offset=getPointOffset(sampleNo,0);
            int res=DataMaths::matrix_inverse(vec, getShape(), offset,
                    temp->getVectorRW(), temp->getShape(), offset, numdpps, h);
            if (res > errorcode) {
                errorcode=res;
#pragma omp critical
                {
                    // I'm not especially concerned which error gets reported
                    // as long as one is
                    errcode=errorcode;
                }
            }
        }
    }
    return errcode;
}

void DataExpanded::setToZero()
{
    // Q: Surely there is a more efficient way to do this????
    // Why is there no memset here? Parallel issues?
    // A: This ensures that memory is touched by the correct thread.
    CHECK_FOR_EX_WRITE;
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    const DataTypes::ValueType::size_type n = getNoValues();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            double* p = &m_data[getPointOffset(sampleNo,dataPointNo)];
            for (int i=0; i<n; ++i)
                p[i] = 0.;
        }
    }
}

void DataExpanded::dump(const std::string fileName) const
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

    // netCDF error handler
    NcError err(NcError::verbose_nonfatal);
    std::string newFileName(esysUtils::appendRankToFileName(fileName,
                                                            mpi_num, mpi_iam));
    NcFile dataFile(newFileName.c_str(), NcFile::Replace);
    if (!dataFile.is_valid())
        throw DataException("DataExpanded::dump: opening of netCDF file for output failed.");
    if (!dataFile.add_att("type_id",2))
        throw DataException("DataExpanded::dump: appending data type to netCDF file failed.");
    if (!dataFile.add_att("rank",rank))
        throw DataException("DataExpanded::dump: appending rank attribute to netCDF file failed.");
    if (!dataFile.add_att("function_space_type",type))
        throw DataException("DataExpanded::dump: appending function space attribute to netCDF file failed.");
    ndims=rank+2;
    if ( rank >0 ) {
        dims[0]=shape[0];
        if (! (ncdims[0] = dataFile.add_dim("d0",shape[0])) )
            throw DataException("DataExpanded::dump: appending ncdim 0 to netCDF file failed.");
    }
    if ( rank >1 ) {
        dims[1]=shape[1];
        if (! (ncdims[1] = dataFile.add_dim("d1",shape[1])) )
            throw DataException("DataExpanded::dump: appending ncdim 1 to netCDF file failed.");
    }
    if ( rank >2 ) {
        dims[2]=shape[2];
        if (! (ncdims[2] = dataFile.add_dim("d2", shape[2])) )
            throw DataException("DataExpanded::dump: appending ncdim 2 to netCDF file failed.");
    }
    if ( rank >3 ) {
        dims[3]=shape[3];
        if (! (ncdims[3] = dataFile.add_dim("d3", shape[3])) )
            throw DataException("DataExpanded::dump: appending ncdim 3 to netCDF file failed.");
    }
    dims[rank]=getFunctionSpace().getNumDataPointsPerSample();
    if (! (ncdims[rank] = dataFile.add_dim("num_data_points_per_sample", dims[rank])) )
        throw DataException("DataExpanded::dump: appending num_data_points_per_sample to netCDF file failed.");
    dims[rank+1]=getFunctionSpace().getNumSamples();
    if (! (ncdims[rank+1] = dataFile.add_dim("num_samples", dims[rank+1])) )
        throw DataException("DataExpanded::dump: appending num_sample to netCDF file failed.");

    if (getFunctionSpace().getNumSamples() > 0) {
        if (! ( ids = dataFile.add_var("id", ncInt, ncdims[rank+1])) )
            throw DataException("DataExpanded::dump: appending reference id to netCDF file failed.");
        const dim_t* ids_p=getFunctionSpace().borrowSampleReferenceIDs();
        if (! (ids->put(ids_p,dims[rank+1])) )
            throw DataException("DataExpanded::dump: copy reference id  to netCDF buffer failed.");
        if (! ( var = dataFile.add_var("data", ncDouble, ndims, ncdims)) )
            throw DataException("DataExpanded::dump: appending variable to netCDF file failed.");
        if (! (var->put(d_ptr,dims)) )
            throw DataException("DataExpanded::dump: copy data to netCDF buffer failed.");
    }
#else
    throw DataException("DataExpanded::dump: not configured with netCDF. Please contact your installation manager.");
#endif // USE_NETCDF
}

void DataExpanded::setTaggedValue(int tagKey,
                                  const DataTypes::ShapeType& pointshape,
                                  const DataTypes::ValueType& value,
                                  int dataOffset)
{
    CHECK_FOR_EX_WRITE;
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    const DataTypes::ValueType::size_type n = getNoValues();
    const double* in = &value[0+dataOffset];

    if (value.size() != n)
        throw DataException("DataExpanded::setTaggedValue: number of input values does not match number of values per data points.");

#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        if (getFunctionSpace().getTagFromSampleNo(sampleNo) == tagKey) {
            for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
                double* p = &m_data[getPointOffset(sampleNo,dataPointNo)];
                for (int i=0; i<n; ++i)
                    p[i] = in[i];
            }
        }
    }
}


void DataExpanded::reorderByReferenceIDs(dim_t *reference_ids)
{
    CHECK_FOR_EX_WRITE;
    const int numSamples = getNumSamples();
    const DataTypes::ValueType::size_type n = getNoValues() * getNumDPPSample();
    FunctionSpace fs=getFunctionSpace();

    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        const index_t id_in = reference_ids[sampleNo];
        const index_t id = fs.getReferenceIDOfSample(sampleNo);
        if (id != id_in) {
            bool matched=false;
            for (int sampleNo2 = sampleNo+1; sampleNo2 < numSamples; sampleNo2++) {
                if (id == reference_ids[sampleNo2]) {
                    double* p = &m_data[getPointOffset(sampleNo,0)];
                    double* p2 = &m_data[getPointOffset(sampleNo2,0)];
                    for (int i=0; i<n; i++) {
                        const double rtmp=p[i];
                        p[i] = p2[i];
                        p2[i] = rtmp;
                    }
                    reference_ids[sampleNo]=id;
                    reference_ids[sampleNo2]=id_in;
                    matched=true;
                    break;
                }
            }
            if (!matched)
                throw DataException("DataExpanded::reorderByReferenceIDs: unable to reorder sample data by reference ids");
        }
    }
}

DataTypes::ValueType& DataExpanded::getVectorRW()
{
    CHECK_FOR_EX_WRITE;
    return m_data.getData();
}

const DataTypes::ValueType& DataExpanded::getVectorRO() const
{
    return m_data.getData();
}

//void DataExpanded::randomFill(long seed)
//{
//    CHECK_FOR_EX_WRITE;
//
//    DataVector&  dv=getVectorRW();
//    const size_t dvsize=dv.size();
//    esysUtils::randomFillArray(seed, &(dv[0]), dvsize);
//}

}  // end of namespace

