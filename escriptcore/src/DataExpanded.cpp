
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "Data.h"
#include "DataConstant.h"
#include "DataException.h"
#include "DataExpanded.h"
#include "DataVectorOps.h"
#include "DataTagged.h"
#include "Utils.h"

#include <limits>

#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif

using namespace std;
using namespace escript::DataTypes;


#ifdef SLOWSHARECHECK
  #define CHECK_FOR_EX_WRITE do {\
    if (isShared()) {\
        std::ostringstream ss;\
        ss << "Attempt to modify shared object. Line " << __LINE__ << " in "\
           << __FILE__;\
        abort();\
        throw DataException(ss.str());\
    }\
  } while(0);
#else
  #define CHECK_FOR_EX_WRITE
#endif

namespace escript {

DataExpanded::DataExpanded(const WrappedArray& value,
                           const FunctionSpace& what)
  : parent(what,value.getShape())
{
    // initialise the data array for this object
    initialise(what.getNumSamples(),what.getNumDPPSample(), value.isComplex());
    // copy the given value to every data point
    copy(value);
}

DataExpanded::DataExpanded(const DataExpanded& other)
  : parent(other.getFunctionSpace(), other.getShape()),
    m_data_r(other.m_data_r), m_data_c(other.m_data_c)
{
    m_iscompl=other.m_iscompl;
}

DataExpanded::DataExpanded(const DataConstant& other)
  : parent(other.getFunctionSpace(), other.getShape())
{
    // initialise the data array for this object
    initialise(other.getNumSamples(),other.getNumDPPSample(), other.isComplex());
    // DataConstant only has one value, copy this to every data point
    copy(other);
}

DataExpanded::DataExpanded(const DataTagged& other)
  : parent(other.getFunctionSpace(), other.getShape())
{
    // initialise the data array for this object
    initialise(other.getNumSamples(),other.getNumDPPSample(), other.isComplex());
    // for each data point in this object, extract and copy the corresponding
    // data value from the given DataTagged object
    if (isComplex())
    {
	DataTypes::cplx_t dummy=0;
	#pragma omp parallel for
	for (int i=0; i<m_noSamples; i++) {
	    for (int j=0; j<m_noDataPointsPerSample; j++) {
		try {
		    DataTypes::copyPoint(getTypedVectorRW(dummy), getPointOffset(i,j),
					getNoValues(), other.getTypedVectorRO(dummy),
					other.getPointOffset(i,j));
		} catch (std::exception& e) {
		    cerr << e.what() << endl;
		}
	    }
	}      
      
      
    }
    else
    {
	DataTypes::real_t dummy=0;
	#pragma omp parallel for
	for (int i=0; i<m_noSamples; i++) {
	    for (int j=0; j<m_noDataPointsPerSample; j++) {
		try {
		    DataTypes::copyPoint(getTypedVectorRW(dummy), getPointOffset(i,j),
					getNoValues(), other.getTypedVectorRO(dummy),
					other.getPointOffset(i,j));
		} catch (std::exception& e) {
		    cerr << e.what() << endl;
		}
	    }
	}
    }
}

DataExpanded::DataExpanded(const DataExpanded& other,
                           const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
    // initialise this Data object to the shape of the slice
    initialise(other.getNumSamples(),other.getNumDPPSample(), other.isComplex());
    // copy the data
    DataTypes::RegionLoopRangeType region_loop_range =
                                    DataTypes::getSliceRegionLoopRange(region);
    if (isComplex())
    {
	DataTypes::cplx_t dummy=0;
	#pragma omp parallel for
	for (int i=0; i<m_noSamples; i++) {
	    for (int j=0; j<m_noDataPointsPerSample; j++) {
		try {
		    DataTypes::copySlice(getTypedVectorRW(dummy), getShape(),
					getPointOffset(i,j), other.getTypedVectorRO(dummy),
					other.getShape(),
					other.getPointOffset(i,j),
					region_loop_range);
		} catch (std::exception& e) {
		    cerr << e.what() << endl;
		}
	    }
	}      
    }
    else
    {
	DataTypes::real_t dummy=0;
	#pragma omp parallel for
	for (int i=0; i<m_noSamples; i++) {
	    for (int j=0; j<m_noDataPointsPerSample; j++) {
		try {
		    DataTypes::copySlice(getTypedVectorRW(dummy), getShape(),
					getPointOffset(i,j), other.getTypedVectorRO(dummy),
					other.getShape(),
					other.getPointOffset(i,j),
					region_loop_range);
		} catch (std::exception& e) {
		    cerr << e.what() << endl;
		}
	    }
	}
    }
}

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::RealVectorType &data)
  : parent(what,shape)
{
    ESYS_ASSERT(data.size()%getNoValues()==0,
                 "DataExpanded Constructor - size of supplied data is not a multiple of shape size.");

    if (data.size() == getNoValues()) {
        RealVectorType& vec=m_data_r;
        // create the view of the data
        initialise(what.getNumSamples(),what.getNumDPPSample(), false);
        // now we copy this value to all elements
        for (int i=0; i<getLength();) {
            for (unsigned int j=0;j<getNoValues();++j,++i) {
                vec[i]=data[j];
            }
        }
    } else {
        // copy the data in the correct format
        m_data_r = data;
    }
}

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::CplxVectorType &data)
  : parent(what,shape)
{
    ESYS_ASSERT(data.size()%getNoValues()==0,
                 "DataExpanded Constructor - size of supplied data is not a multiple of shape size.");

    if (data.size() == getNoValues()) {
        CplxVectorType& vec=m_data_c;
        // create the view of the data
        initialise(what.getNumSamples(),what.getNumDPPSample(), true);
        // now we copy this value to all elements
        for (int i=0; i<getLength();) {
            for (unsigned int j=0;j<getNoValues();++j,++i) {
                vec[i]=data[j];
            }
        }
    } else {
        // copy the data in the correct format
        m_data_c = data;
    }
}


DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::real_t v)
  : parent(what,shape)
{
    initialise(what.getNumSamples(),what.getNumDPPSample(), false);
    DataTypes::RealVectorType& vec=m_data_r;
    // now we copy this value to all elements
    const int L=getLength();
#pragma omp parallel for
    for (int i=0; i<L; ++i) {
        vec[i]=v;
    }
}

DataExpanded::DataExpanded(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::cplx_t v)
  : parent(what,shape)
{
    initialise(what.getNumSamples(),what.getNumDPPSample(), true);
    DataTypes::CplxVectorType& vec=m_data_c;
    
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

DataAbstract* DataExpanded::deepCopy() const
{
    return new DataExpanded(*this);
}

DataAbstract*
DataExpanded::zeroedCopy() const
{
    DataExpanded* p=0;
    if (isComplex())
    {
        p=new DataExpanded(this->getFunctionSpace(), this->getShape(), DataTypes::cplx_t(0));
    }
    else
    {
        p=new DataExpanded(this->getFunctionSpace(), this->getShape(), DataTypes::real_t(0));
    }   
  return p;
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

    if (value->isComplex()!=isComplex())
    {
	throw DataException("Programmer Error: object and new value must be the same complexity.");
    }

    if (isComplex())
    {
	// copy the data from the slice into this object
	DataTypes::CplxVectorType& vec=getVectorRWC();
	const ShapeType& mshape=getShape();
	const DataTypes::CplxVectorType& tVec=tempDataExp->getVectorROC();
	const ShapeType& tShape=tempDataExp->getShape();
    #pragma omp parallel for
	for (int i=0; i<m_noSamples; i++) {
	    for (int j=0; j<m_noDataPointsPerSample; j++) {
		DataTypes::copySliceFrom(vec, mshape, getPointOffset(i,j), tVec,
					tShape, tempDataExp->getPointOffset(i,j),
					region_loop_range);
	    }
	}      
    }
    else
    {
	// copy the data from the slice into this object
	DataTypes::RealVectorType& vec=getVectorRW();
	const ShapeType& mshape=getShape();
	const DataTypes::RealVectorType& tVec=tempDataExp->getVectorRO();
	const ShapeType& tShape=tempDataExp->getShape();
    #pragma omp parallel for
	for (int i=0; i<m_noSamples; i++) {
	    for (int j=0; j<m_noDataPointsPerSample; j++) {
		DataTypes::copySliceFrom(vec, mshape, getPointOffset(i,j), tVec,
					tShape, tempDataExp->getPointOffset(i,j),
					region_loop_range);
	    }
	}
    }
}

void DataExpanded::copy(const DataConstant& value)
{
    ESYS_ASSERT(checkShape(getShape(), value.getShape()),
                 createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.", value.getShape(), getShape()));
    if (isComplex())
    {
	if (value.isComplex())
	{
	    // copy a single value to every data point in this object
	    #pragma omp parallel for
	    for (int i=0; i<m_noSamples; i++) {
		for (int j=0; j<m_noDataPointsPerSample; j++) {
		    DataTypes::copyPoint(getTypedVectorRW((cplx_t)(0)), getPointOffset(i,j),
					getNoValues(), value.getTypedVectorRO((cplx_t)(0)), 0);
		}
	    }	    
	}
	else	// value is real
	{
	    throw DataException("Programming error - DataExpanded::copy source and target must be the same complexity.");	  
	}
    }
    else
    {
	if (value.isComplex())
	{
	    throw DataException("Programming error - DataExpanded::copy source and target must be the same complexity.");	  	  
	}
	else
	{
	    real_t dummy=0;
	    // copy a single value to every data point in this object
	    #pragma omp parallel for
	    for (int i=0; i<m_noSamples; i++) {
		for (int j=0; j<m_noDataPointsPerSample; j++) {
		    DataTypes::copyPoint(getTypedVectorRW(dummy), getPointOffset(i,j),
					getNoValues(), value.getTypedVectorRO(dummy), 0);
		}
	    }
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

void DataExpanded::initialise(int noSamples, int noDataPointsPerSample, bool cplx)
{
    this->m_iscompl=cplx;
    if (noSamples==0) //retain the default empty object
        return;

    if (cplx)
    {
	// resize data array to the required size
	m_data_c.resize(noSamples*noDataPointsPerSample*getNoValues(), 0.0, noDataPointsPerSample*getNoValues());      
    }
    else
    {
	// resize data array to the required size
	m_data_r.resize(noSamples*noDataPointsPerSample*getNoValues(), 0.0, noDataPointsPerSample*getNoValues());
    }
}

bool
DataExpanded::hasNaN() const
{
  bool haveNaN=false;
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
	  if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag()))
	  {
	      #pragma omp critical 
	      {
		  haveNaN=true;
	      }
	  }
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
	  if (std::isnan(m_data_r[i]))
	  {
	      #pragma omp critical 
	      {
		  haveNaN=true;
	      }
	  }
      }
  }
  return haveNaN;
}


void
DataExpanded::replaceNaN(DataTypes::real_t value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
	if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag()))  
	{
	  m_data_c[i] = value;
	}
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
	if (std::isnan(m_data_r[i]))  
	{
	  m_data_r[i] = value;
	}
      }    
  }
}

void
DataExpanded::replaceNaN(DataTypes::cplx_t value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
	if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag())) 
	{
	  m_data_c[i] = value;
	}
      }
  }
  else
  {
      complicate();
      replaceNaN(value);
  }
}


bool
DataExpanded::hasInf() const
{
  bool haveNaN=false;
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
	  if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag()))
	  {
	      #pragma omp critical 
	      {
		  haveNaN=true;
	      }
	  }
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
	  if (std::isinf(m_data_r[i]))
	  {
	      #pragma omp critical 
	      {
		  haveNaN=true;
	      }
	  }
      }
  }
  return haveNaN;
}


void
DataExpanded::replaceInf(DataTypes::real_t value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
	if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag()))  
	{
	  m_data_c[i] = value;
	}
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
	if (std::isinf(m_data_r[i]))  
	{
	  m_data_r[i] = value;
	}
      }    
  }
}

void
DataExpanded::replaceInf(DataTypes::cplx_t value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
	if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag())) 
	{
	  m_data_c[i] = value;
	}
      }
  }
  else
  {
      complicate();
      replaceInf(value);
  }
}



string DataExpanded::toString() const
{
    stringstream ss;
    FunctionSpace fs=getFunctionSpace();

    int offset=0;
    for (int i=0; i<m_noSamples; i++) {
        for (int j=0; j<m_noDataPointsPerSample; j++) {
            offset = getPointOffset(i,j);
            stringstream suffix;
            suffix << "( id: " << i << ", ref: "
                   << fs.getReferenceIDOfSample(i) << ", pnt: " << j << ")";
            ss << (isComplex()?
		    DataTypes::pointToString(getTypedVectorRO((cplx_t)0), getShape(), offset, suffix.str())
		   :
		    DataTypes::pointToString(getTypedVectorRO((real_t)0), getShape(), offset, suffix.str()));
            if (!(i==(m_noSamples-1) && j==(m_noDataPointsPerSample-1))) {
                ss << endl;
            }
        }
    }
    string result(ss.str());
    if (result.empty())
        return "(data contains no samples)\n";

    return result;
}

DataTypes::RealVectorType::size_type DataExpanded::getPointOffset(int sampleNo,
                                                        int dataPointNo) const
{
    DataTypes::RealVectorType::size_type blockSize=getNoValues();
    ESYS_ASSERT((isComplex()?
		  ((sampleNo >= 0) && (dataPointNo >= 0) && (m_data_c.size() > 0))
		:
		  ((sampleNo >= 0) && (dataPointNo >= 0) && (m_data_r.size() > 0))), 
	       "(DataBlocks2D) Index value out of range.");
    DataTypes::RealVectorType::size_type temp=(sampleNo*m_noDataPointsPerSample+dataPointNo)*blockSize;
    ESYS_ASSERT((isComplex()?
		  (temp <= (m_data_c.size()-blockSize))
		:
		  (temp <= (m_data_r.size()-blockSize))), "Index value out of range.");

    return temp;
}


void DataExpanded::complicate()
{
    if (!isComplex())
    {
        fillComplexFromReal(m_data_r, m_data_c);
        this->m_iscompl=true;
        m_data_r.resize(0,0,1);
    }
}


DataTypes::RealVectorType::size_type DataExpanded::getLength() const
{
    return std::max(m_data_c.size(), m_data_r.size());  
}

void DataExpanded::copyToDataPoint(int sampleNo, int dataPointNo, const DataTypes::cplx_t value)
{
    if (!isComplex())
    {
	throw DataException("Programming error - attempt to set complex value on real data.");
    }
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

        DataTypes::CplxVectorType::size_type offset = getPointOffset(sampleNo, dataPointNo);
        DataTypes::CplxVectorType& vec = getTypedVectorRW(cplx_t(0));
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


void DataExpanded::copyToDataPoint(int sampleNo, int dataPointNo, const DataTypes::real_t value)
{
    if (isComplex())
    {
	copyToDataPoint(sampleNo, dataPointNo, cplx_t(value));
	return;
    }
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

        DataTypes::RealVectorType::size_type offset = getPointOffset(sampleNo, dataPointNo);
        DataTypes::RealVectorType& vec = getVectorRW();
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

	if (isComplex())
	{
	    DataTypes::CplxVectorType::size_type offset = getPointOffset(sampleNo, dataPointNo);
	    DataTypes::CplxVectorType& vec = getTypedVectorRW(cplx_t(0));
	    vec.copyFromArrayToOffset(value, offset, 1);
	}
	else
	{
	    DataTypes::RealVectorType::size_type offset = getPointOffset(sampleNo, dataPointNo);
	    DataTypes::RealVectorType& vec = getTypedVectorRW(real_t(0));
	    vec.copyFromArrayToOffset(value, offset, 1);
	}
    }
}

void DataExpanded::symmetric(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev = dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::symmetric: casting to DataExpanded failed (probably a programming error).");

    const ShapeType& shape = getShape();
    const ShapeType& evShape = temp_ev->getShape();
    if (isComplex())
    {
	const DataTypes::CplxVectorType& vec = getTypedVectorRO((DataTypes::cplx_t)0);
	DataTypes::CplxVectorType& evVec = temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::symmetric(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo));
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType& vec = getTypedVectorRO(0.0);
	DataTypes::RealVectorType& evVec = temp_ev->getTypedVectorRW(0.0);
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::symmetric(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo));
	    }
	}
    }
}

void DataExpanded::antisymmetric(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::antisymmetric: casting to DataExpanded failed (probably a programming error).");

    const ShapeType& shape = getShape();
    const ShapeType& evShape = temp_ev->getShape();
    if (isComplex())
    {
	const DataTypes::CplxVectorType& vec = getTypedVectorRO((DataTypes::cplx_t)0);
	DataTypes::CplxVectorType& evVec = temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) 
	{
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++)
	    {
		escript::antisymmetric(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo));
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType& vec = getTypedVectorRO(0.0);
	DataTypes::RealVectorType& evVec = temp_ev->getTypedVectorRW(0.0);
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++)
	    {
		escript::antisymmetric(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo));
	    }
	}
    }    
}


void DataExpanded::hermitian(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev = dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::hermitian: casting to DataExpanded failed (probably a programming error).");
    if (!isComplex() || !temp_ev->isComplex())
    {
	throw DataException("DataExpanded::hermitian: do not call this method with real data");
    }
    const ShapeType& shape = getShape();
    const ShapeType& evShape = temp_ev->getShape();
    const DataTypes::CplxVectorType& vec = getTypedVectorRO((DataTypes::cplx_t)0);
    DataTypes::CplxVectorType& evVec = temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
    #pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
	    escript::hermitian(vec, shape,
		    getPointOffset(sampleNo,dataPointNo), evVec, evShape,
		    ev->getPointOffset(sampleNo,dataPointNo));
	}
    }
}

void DataExpanded::antihermitian(DataAbstract* ev)
{
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    DataExpanded* temp_ev=dynamic_cast<DataExpanded*>(ev);
    if (!temp_ev)
        throw DataException("DataExpanded::antihermitian: casting to DataExpanded failed (probably a programming error).");
    if (!isComplex() || !temp_ev->isComplex())
    {
	throw DataException("DataExpanded::antihermitian: do not call this method with real data");
    }
    const ShapeType& shape = getShape();
    const ShapeType& evShape = temp_ev->getShape();
    const DataTypes::CplxVectorType& vec = getTypedVectorRO((DataTypes::cplx_t)0);
    DataTypes::CplxVectorType& evVec = temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
    #pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) 
    {
	for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++)
	{
	    escript::antihermitian(vec, shape,
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
    const ShapeType& shape=getShape();
    const ShapeType& evShape=temp_ev->getShape(); 
    if (isComplex())
    {
	const DataTypes::CplxVectorType& vec=getVectorROC();
	DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();

    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::trace(vec, shape, getPointOffset(sampleNo,dataPointNo),
			evVec, evShape, ev->getPointOffset(sampleNo,dataPointNo),
			axis_offset);
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType& vec=getVectorRO();
	DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();

    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::trace(vec, shape, getPointOffset(sampleNo,dataPointNo),
			evVec, evShape, ev->getPointOffset(sampleNo,dataPointNo),
			axis_offset);
	    }
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
    const ShapeType& shape = getShape();
    if (isComplex())
    {
	const DataTypes::CplxVectorType& vec = getVectorROC();
	DataTypes::CplxVectorType& evVec = temp_ev->getVectorRWC();
	const ShapeType& evShape = temp_ev->getShape();
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::transpose(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo), axis_offset);
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType& vec = getVectorRO();
	DataTypes::RealVectorType& evVec = temp_ev->getVectorRW();
	const ShapeType& evShape = temp_ev->getShape();
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::transpose(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo), axis_offset);
	    }
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
    const ShapeType& shape=getShape();    
    const ShapeType& evShape=temp_ev->getShape();
    if (isComplex())
    {
	const DataTypes::CplxVectorType& vec=getVectorROC();
	DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();

    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::swapaxes(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo), axis0, axis1);
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType& vec=getVectorRO();
	DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();

    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::swapaxes(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo), axis0, axis1);
	    }
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
    const ShapeType& evShape=temp_ev->getShape();
    const ShapeType& shape=getShape();
    if (isComplex())
    {
	const DataTypes::CplxVectorType& vec=getVectorROC();
	DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();

    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::eigenvalues(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo));
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType& vec=getVectorRO();
	DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();

    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		escript::eigenvalues(vec, shape,
			getPointOffset(sampleNo,dataPointNo), evVec, evShape,
			ev->getPointOffset(sampleNo,dataPointNo));
	    }
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

    const DataTypes::RealVectorType& vec = getVectorRO();
    const ShapeType& shape = getShape();
    DataTypes::RealVectorType& evVec = temp_ev->getVectorRW();
    const ShapeType& evShape = temp_ev->getShape();
    DataTypes::RealVectorType& VVec = temp_V->getVectorRW();
    const ShapeType& VShape = temp_V->getShape();
#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
            escript::eigenvalues_and_eigenvectors(vec, shape,
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
    const DataTypes::RealVectorType& vec=m_data_r;
    int errcode=0;
#pragma omp parallel
    {
        int errorcode=0;
        LapackInverseHelper h(getShape()[0]);
#pragma omp for
        for (int sampleNo = 0; sampleNo < numSamples; sampleNo++)
        {
            // not sure I like all those virtual calls to getPointOffset
            DataTypes::RealVectorType::size_type offset=getPointOffset(sampleNo,0);
            int res=escript::matrix_inverse(vec, getShape(), offset,
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
    if (isComplex())
    {
	const DataTypes::CplxVectorType::size_type n = getNoValues();
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		DataTypes::cplx_t* p = &m_data_c[getPointOffset(sampleNo,dataPointNo)];
		for (int i=0; i<n; ++i)
		    p[i] = 0.;
	    }
	}
    }
    else
    {
	const DataTypes::RealVectorType::size_type n = getNoValues();
    #pragma omp parallel for
	for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
	    for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
		double* p = &m_data_r[getPointOffset(sampleNo,dataPointNo)];
		for (int i=0; i<n; ++i)
		    p[i] = 0.;
	    }
	}
    }
}

#ifdef ESYS_HAVE_HDF5
void DataExpanded::dump_hdf5(const H5::Group h5_grp) const
{
    uint rank = getRank();
    int fs_type=  getFunctionSpace().getTypeCode();
    const DataTypes::ShapeType& shape = getShape();
    const dim_t* ids_p=getFunctionSpace().borrowSampleReferenceIDs();
    if (isComplex())
    {
        throw DataException("Error - DataExpanded::dump_hdf5: complex data are not supported. Split into real and imaginary part.");
    }
#ifdef ESYS_MPI
   /* Serialize I/O */
   JMPI mpiInfo = getFunctionSpace().getDomain()->getMPI();
   const int mpi_iam = mpiInfo->rank;
   const int mpi_num = mpiInfo->size;
   MPI_Comm comm = mpiInfo->comm;
   MPI_Status status;
   int dummy = 0;
   if (mpi_iam > 0)
       MPI_Recv(&dummy, 1, MPI_INT, mpi_iam-1, 81802, comm, &status);
#endif

    try
    {
        // .... add meta data ............
        uint h5_shape[DataTypes::maxRank]; // dataset dimensions
        for (uint i = 0; i < rank; i++) {
            h5_shape[i]= shape[i];
        }
        hsize_t h5_shape_dims[1] = {static_cast<hsize_t>(rank)};
        H5::DataSet h5_dmeta = h5_grp.createDataSet("Meta", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_shape_dims ) );
        h5_dmeta.write( h5_shape, H5::PredType::NATIVE_UINT);
        // data type
        hsize_t h5_typeid_dims[1] = { 1 };
        int h5_type_id[1] = { 2 };
        H5::Attribute h5_typeid_attr = h5_dmeta.createAttribute("type_id", H5::PredType::NATIVE_INT, H5::DataSpace(1,h5_typeid_dims ) );
        h5_typeid_attr.write( H5::PredType::NATIVE_INT , h5_type_id );

        uint h5_rank[1] = { static_cast<uint>(rank) };
        H5::Attribute h5_rank_attr = h5_dmeta.createAttribute("rank", H5::PredType::NATIVE_UINT, H5::DataSpace(1,h5_typeid_dims ) );
        h5_rank_attr.write( H5::PredType::NATIVE_UINT , h5_rank );

        int dfs_type[1] = { fs_type };
        H5::Attribute h5_fs_type_attr = h5_dmeta.createAttribute("function_space_type", H5::PredType::NATIVE_INT, H5::DataSpace(1,h5_typeid_dims) );
        h5_fs_type_attr.write( H5::PredType::NATIVE_INT , dfs_type );
        // .... end meta data ............
        // ... add index of samples ....
        hsize_t h5_id_dims[1] = { static_cast<hsize_t>(getFunctionSpace().getNumSamples()) };
        H5::DataSpace h5_dataspace_ids(1 , h5_id_dims );
        #ifdef ESYS_INDEXTYPE_LONG
            H5::DataSet h5_dataset_ids = h5_grp.createDataSet("sample_id", H5::PredType::NATIVE_LONG, h5_dataspace_ids );
            h5_dataset_ids.write(ids_p, H5::PredType::NATIVE_LONG);
        #else
            H5::DataSet h5_dataset_ids = h5_grp.createDataSet("sample_id", H5::PredType::NATIVE_INT, h5_dataspace_ids );
            h5_dataset_ids.write(ids_p, H5::PredType::NATIVE_INT);
        #endif
        // ... add data ....
        const double* d_ptr=&(m_data_r[0]);
        hsize_t h5_data_length[1] = {  static_cast<hsize_t>(m_data_r.size()) };
        H5::DataSet h5_dataset_data = h5_grp.createDataSet("data", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1 , h5_data_length ) );
        h5_dataset_data.write(d_ptr, H5::PredType::NATIVE_DOUBLE);
    }
    // catch failure caused by the H5File operations
    catch (H5::Exception& error)
    {
        #ifdef ESYS_MPI
           if ( mpi_iam < mpi_num-1 ) MPI_Send(&dummy, 1, MPI_INT, mpi_iam+1, 81802, comm);
        #endif
        error.printErrorStack();
        throw DataException("Error - DataExpanded:: creating HDF5 file failed.");
    }
    #ifdef ESYS_MPI
       if ( mpi_iam < mpi_num-1 ) MPI_Send(&dummy, 1, MPI_INT, mpi_iam+1, 81802, comm);
    #endif
}
#endif

void DataExpanded::setTaggedValue(int tagKey,
                                  const DataTypes::ShapeType& pointshape,
                                  const DataTypes::RealVectorType& value,
                                  int dataOffset)
{
    CHECK_FOR_EX_WRITE;
    if (isComplex())
    {
        CplxVectorType tv;
	fillComplexFromReal(value, tv);
	setTaggedValue(tagKey, pointshape, tv, dataOffset);
        return;
    }
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    const DataTypes::RealVectorType::size_type n = getNoValues();
    const real_t* in = &value[0+dataOffset];

    if (value.size() != n)
        throw DataException("DataExpanded::setTaggedValue: number of input values does not match number of values per data points.");

#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        if (getFunctionSpace().getTagFromSampleNo(sampleNo) == tagKey) {
            for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
                real_t* p = &m_data_r[getPointOffset(sampleNo,dataPointNo)];
                for (int i=0; i<n; ++i)
                    p[i] = in[i];
            }
        }
    }
}


void DataExpanded::setTaggedValue(int tagKey,
                                  const DataTypes::ShapeType& pointshape,
                                  const DataTypes::CplxVectorType& value,
                                  int dataOffset)
{
    CHECK_FOR_EX_WRITE;
    if (!isComplex())
    {
	throw DataException("Programming Error - Attempt to set a complex value on a real object.");
    }
    const int numSamples = getNumSamples();
    const int numDataPointsPerSample = getNumDPPSample();
    const DataTypes::CplxVectorType::size_type n = getNoValues();
    const DataTypes::cplx_t* in = &value[0+dataOffset];

    if (value.size() != n)
        throw DataException("DataExpanded::setTaggedValue: number of input values does not match number of values per data points.");

#pragma omp parallel for
    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        if (getFunctionSpace().getTagFromSampleNo(sampleNo) == tagKey) {
            for (int dataPointNo = 0; dataPointNo < numDataPointsPerSample; dataPointNo++) {
                cplx_t* p = &m_data_c[getPointOffset(sampleNo,dataPointNo)];
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
    const DataTypes::RealVectorType::size_type n = getNoValues() * getNumDPPSample();
    FunctionSpace fs=getFunctionSpace();

    for (int sampleNo = 0; sampleNo < numSamples; sampleNo++) {
        const index_t id_in = reference_ids[sampleNo];
        const index_t id = fs.getReferenceIDOfSample(sampleNo);
        if (id != id_in) {
            bool matched=false;
            for (int sampleNo2 = sampleNo+1; sampleNo2 < numSamples; sampleNo2++) {
                if (id == reference_ids[sampleNo2]) {
                    double* p = &m_data_r[getPointOffset(sampleNo,0)];
                    double* p2 = &m_data_r[getPointOffset(sampleNo2,0)];
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

DataTypes::RealVectorType& DataExpanded::getVectorRW()
{
    CHECK_FOR_EX_WRITE;
    return m_data_r;
}

const DataTypes::RealVectorType& DataExpanded::getVectorRO() const
{
    return m_data_r;
}

DataTypes::CplxVectorType& DataExpanded::getVectorRWC()
{
    CHECK_FOR_EX_WRITE;
    return m_data_c;
}

const DataTypes::CplxVectorType& DataExpanded::getVectorROC() const
{
    return m_data_c;
}

DataTypes::RealVectorType& DataExpanded::getTypedVectorRW(DataTypes::real_t dummypar)
{
    CHECK_FOR_EX_WRITE;
    return m_data_r;
}

const DataTypes::RealVectorType& DataExpanded::getTypedVectorRO(DataTypes::real_t dummypar) const
{
    return m_data_r;
}

DataTypes::CplxVectorType& DataExpanded::getTypedVectorRW(DataTypes::cplx_t dummypar)
{
    CHECK_FOR_EX_WRITE;
    return m_data_c;
}

const DataTypes::CplxVectorType& DataExpanded::getTypedVectorRO(DataTypes::cplx_t dummypar) const
{
    return m_data_c;
}


//void DataExpanded::randomFill(long seed)
//{
//    CHECK_FOR_EX_WRITE;
//
//    DataVector&  dv=getVectorRW();
//    const size_t dvsize=dv.size();
//    randomFillArray(seed, &(dv[0]), dvsize);
//}

}  // end of namespace

