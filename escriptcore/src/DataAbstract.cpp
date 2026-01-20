
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "DataAbstract.h"
#include "Data.h" // So we can update the shared status when things change
#include "DataException.h"
#include "DataLazy.h"


using namespace std;

namespace escript {

// The boost methods enable_shared_from_this return a shared_ptr managing the object
// when all you have is the object. It explicitly fails in the case where
// you haven't made a shared_ptr for this object yet.
// _Currently_ we need this behaviour, hence the exception squashing.
// An exception would be thrown in two circumstances:
//    1. The object doesn't have a shared_ptr attached yet.
//    2. All shared_ptrs have let go and the object is in the process of being destroyed.
// An attempt to getPtr() in the second case is doomed anyway.
// 
// Use of something equivalent to boost_1_39's make_shared used elsewhere might remove
// the need for the hack.
DataAbstract_ptr DataAbstract::getPtr()
{
  try
  {
      return shared_from_this();
  }
  catch (boost::bad_weak_ptr& p)     
  {
      return DataAbstract_ptr(this);
  }
}

const_DataAbstract_ptr DataAbstract::getPtr() const 
{
  try
  {
      return shared_from_this();
  }
  catch (boost::bad_weak_ptr& p)     
  {
      return const_DataAbstract_ptr(this);
  }
}

bool
DataAbstract::isLazy() const
{
    return (dynamic_cast<const DataLazy*>(this)!=0);
}

bool
DataAbstract::isComplex() const
{
    return m_iscompl;
}


DataAbstract::DataAbstract(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty, bool isCplx):
    m_noSamples(what.getNumSamples()),
    m_noDataPointsPerSample(what.getNumDPPSample()),
    m_iscompl(isCplx),
    m_functionSpace(what),
    m_shape(shape),
    m_novalues(DataTypes::noValues(shape)),
    m_rank(DataTypes::getRank(shape))
{
#ifdef EXWRITECHK
    exclusivewritecalled=false;
#endif  
  
    m_isempty=isDataEmpty;
    if (m_rank>ESCRIPT_MAX_DATA_RANK)
    {
    ostringstream os;
        os << "Error - Attempt to create a rank " << m_rank 
       << " object. The maximum rank is " << ESCRIPT_MAX_DATA_RANK << ".";
     throw DataException(os.str());
    }
}

DataAbstract::~DataAbstract() 
{
}


void
DataAbstract::operandCheck(const DataAbstract& right) const
{
    if ((right.getNumDPPSample()!=getNumDPPSample()) ||
    (right.getNumSamples()!=getNumSamples()) ||
    (right.getFunctionSpace()!=getFunctionSpace())) {
      stringstream temp;
      temp << "Error - Right hand argument sample shape or function space "
       << "incompatible with left." << endl 
       << "LHS: (" << getNumSamples() << ","
       << getNumDPPSample() << ") " << getFunctionSpace().toString() 
       << endl
       << "RHS: (" << right.getNumSamples() << "," 
       << right.getNumDPPSample() << ") " 
       << right.getFunctionSpace().toString();
      throw DataException(temp.str());
    }

    //
    // Check the shape of the point data, a rank of 0(scalar) is okay
    if (!((right.getRank()==0) || (getRank()==0) || 
      (right.getShape()==getShape())))
      {
        stringstream temp;
    temp << "Error - Right hand argument point data shape: " 
         << DataTypes::shapeToString(right.getShape())
         << " doesn't match left: " 
         << DataTypes::shapeToString(getShape());
    throw DataException(temp.str());
      }
}

#ifdef ESYS_HAVE_HDF5
void
DataAbstract::dump_hdf5(const H5::Group h5_grp) const
{
    throw DataException("Error - DataAbstract::dump_hdf5: not implemented.");
}
#endif

DataTypes::real_t*
DataAbstract::getSampleDataByTag(int tag, DataTypes::real_t dummy)
{
    throw DataException("Error - DataAbstract::getSampleDataByTag: Data type does not have tag values.");
}


DataTypes::cplx_t*
DataAbstract::getSampleDataByTag(int tag, DataTypes::cplx_t dummy)
{
    throw DataException("Error - DataAbstract::getSampleDataByTag_C: Data type does not have complex tag values.");
}

size_t
DataAbstract::getTagCount() const
{
    return 0;
}



void  
DataAbstract::setTaggedValue(int tagKey,
           const DataTypes::ShapeType& pointshape,
           const DataTypes::RealVectorType& value,
           int dataOffset)
{
    throw DataException("Error - DataAbstract::setTaggedValue: Data type does not have tag values.");
}

void  
DataAbstract::setTaggedValue(int tagKey,
           const DataTypes::ShapeType& pointshape,
           const DataTypes::CplxVectorType& value,
           int dataOffset)
{
    throw DataException("Error - DataAbstract::setTaggedValue: Data type does not have tag values.");
}

int
DataAbstract::getTagNumber(int dpno)
{
    throw DataException("Error - DataAbstract::getTagNumber: Data type cannot be accessed by tag values.");
    return (0);
}

void
DataAbstract::copyToDataPoint(const int sampleNo, const int dataPointNo, const DataTypes::real_t value)
{
    throw DataException("Error - DataAbstract::copying data from double value to a single data point is not supported.");
}

void
DataAbstract::copyToDataPoint(const int sampleNo, const int dataPointNo, const DataTypes::cplx_t value)
{
    throw DataException("Error - DataAbstract::copying data from double value to a single data point is not supported.");
}


void
DataAbstract::copyToDataPoint(const int sampleNo, const int dataPointNo, const WrappedArray& value)
{
    throw DataException("Error - DataAbstract::copying data from WrappedArray objects to a single data point is not supported.");
}


void
DataAbstract::symmetric(DataAbstract* ev) 
{
    throw DataException("Error - DataAbstract::symmetric is not supported.");
}

void
DataAbstract::antisymmetric(DataAbstract* ev) 
{
    throw DataException("Error - DataAbstract::antisymmetric is not supported.");
}

void
DataAbstract::hermitian(DataAbstract* ev) 
{
    throw DataException("Error - DataAbstract::hermitian is not supported.");
}

void
DataAbstract::antihermitian(DataAbstract* ev) 
{
    throw DataException("Error - DataAbstract::antihermitian is not supported.");
}

void
DataAbstract::trace(DataAbstract* ev, int axis_offset) 
{
    throw DataException("Error - DataAbstract::trace is not supported.");
}

void
DataAbstract::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
    throw DataException("Error - DataAbstract::component swapaxes is not supported.");
}
void
DataAbstract::transpose(DataAbstract* ev, int axis_offset) 
{
    throw DataException("Error - DataAbstract::transpose is not supported.");
}

void
DataAbstract::eigenvalues(DataAbstract* ev) 
{
    throw DataException("Error - DataAbstract::eigenvalues is not supported.");
}
void
DataAbstract::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
    throw DataException("Error - DataAbstract::eigenvalues_and_eigenvectors is not supported.");

}

int
DataAbstract::matrixInverse(DataAbstract* out) const
{
   throw DataException("Error - DataAbstract::matrixInverse is not supported.");
}

void
DataAbstract::setToZero() 
{
    throw DataException("Error - DataAbstract:: cannot set values to zero.");
}

void
DataAbstract::reorderByReferenceIDs(DataTypes::dim_t* reference_ids)
{
    throw DataException("Error - DataAbstract:: cannot reorder by reference ids.");
}

void DataAbstract::complicate()
{
    throw DataException("This type does not support converting to complex.");
}


}  // end of namespace

