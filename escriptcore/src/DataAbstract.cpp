
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
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


#include "DataAbstract.h"
#include "DataException.h"
#include "DataLazy.h"

#include "Data.h" // So we can update the shared status when things change

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
  catch (boost::bad_weak_ptr p)     
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
  catch (boost::bad_weak_ptr p)     
  {
      return const_DataAbstract_ptr(this);
  }
}


// Warning - this method uses .use_count() which the boost doco labels inefficient.
// If this method needs to be called in debug contexts, we may need to do some
// timing experiments to determine how inefficient and possibly switch over to
// invasive pointers which can answer these questions faster
bool DataAbstract::checkNoSharing() const
{

  return !m_lazyshared && (m_owners.size()<2);

/*  if (_internal_weak_this.expired())  // there is no shared_ptr for this object yet
  {
    return true;
  }
  if (shared_from_this().use_count()==2)    // shared_from_this will increase the ref count
  {                     // which is the reason .unique is no use.
    return true;
  }
std::cerr << "-<"<<shared_from_this().use_count() << ">-" << endl;
  return false;*/
}

bool
DataAbstract::isLazy() const
{
    return (dynamic_cast<const DataLazy*>(this)!=0);
}



DataAbstract::DataAbstract(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty):
    m_lazyshared(false),
    m_noSamples(what.getNumSamples()),
    m_noDataPointsPerSample(what.getNumDPPSample()),
    m_functionSpace(what),
    m_shape(shape),
    m_novalues(DataTypes::noValues(shape)),
    m_rank(DataTypes::getRank(shape))

{
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

void
DataAbstract::dump(const std::string fileName) const
{
    throw DataException("Error - DataAbstract::dump: not implemented.");
}



DataAbstract::ValueType::value_type*
DataAbstract::getSampleDataByTag(int tag)
{
    throw DataException("Error - DataAbstract::getSampleDataByTag: Data type does not have tag values.");
}

size_t
DataAbstract::getTagCount() const
{
    return 0;
}



void  
DataAbstract::setTaggedValue(int tagKey,
           const DataTypes::ShapeType& pointshape,
               const DataTypes::ValueType& value,
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
DataAbstract::copyToDataPoint(const int sampleNo, const int dataPointNo, const double value)
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
DataAbstract::nonsymmetric(DataAbstract* ev) 
{
    throw DataException("Error - DataAbstract::nonsymmetric is not supported.");
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
DataAbstract::reorderByReferenceIDs(dim_t *reference_ids)
{
    throw DataException("Error - DataAbstract:: cannot reorder by reference ids.");
}


void DataAbstract::addOwner(Data* d)
{
  for (size_t i=0;i<m_owners.size();++i)
  {
    if (m_owners[i]==d)
    {
        return;
    }
  }
  m_owners.push_back(d);
// cerr << "Adding " << d << " as an owner of " << this << " now O=" << m_owners.size() << endl;
  if (m_owners.size()==2)   // Means it used to be 1 so we need to tell people
  {
    for (size_t i=0;i<m_owners.size();++i)
    {
        m_owners[i]->updateShareStatus(true);
    }
  }
}

void DataAbstract::removeOwner(Data* d)
{
  for (size_t i=0;i<m_owners.size();++i)
  {
    if (m_owners[i]==d)
    {
        m_owners.erase(m_owners.begin()+i,m_owners.begin()+(i+1));  // remove the element
        break;
    }
  }
  if (m_owners.size()==1)   // Means it used to be 2 so we need to tell people
  {
    m_owners[0]->updateShareStatus(isShared());     // could still be lazy shared
  }
}


void DataAbstract::makeLazyShared()
{
    m_lazyshared=true;  // now we need to inform all the owners
    for (size_t i=0;i<m_owners.size();++i)
    {
        m_owners[i]->updateShareStatus(true);
    }
}   


}  // end of namespace
