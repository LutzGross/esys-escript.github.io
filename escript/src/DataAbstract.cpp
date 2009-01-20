
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


#include "DataAbstract.h"
#include "DataException.h"
#include "DataLazy.h"

#include "Data.h"		// So we can update the shared status when things change

using namespace std;

namespace escript {

DataAbstract_ptr DataAbstract::getPtr()
{
  if (_internal_weak_this.expired())
  {
	return DataAbstract_ptr(this);	
  }
  else
  {
	return shared_from_this();
  }
}

const_DataAbstract_ptr DataAbstract::getPtr() const 
{
  if (_internal_weak_this.expired())
  {
	return const_DataAbstract_ptr(this);
  }
  else
  {
	return shared_from_this();
  }
}


// Warning - this method uses .use_count() which the boost doco labels inefficient.
// If this method needs to be called in debug contexts, we may need to do some
// timing experiments to determine how inefficient and possibly switch over to
// invasive pointers which can answer these questions faster
bool DataAbstract::checkNoSharing() const
{

  return !m_lazyshared && (m_owners.size()<2);

/*  if (_internal_weak_this.expired())	// there is no shared_ptr for this object yet
  {
	return true;
  }
  if (shared_from_this().use_count()==2)	// shared_from_this will increase the ref count
  {						// which is the reason .unique is no use.
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
    m_noSamples(what.getNumSamples()),
    m_noDataPointsPerSample(what.getNumDPPSample()),
    m_functionSpace(what),
    m_shape(shape),
    m_novalues(DataTypes::noValues(shape)),
    m_rank(DataTypes::getRank(shape)),
    m_lazyshared(false)

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
void
DataAbstract::setToZero() 
{
    throw DataException("Error - DataAbstract:: cannot set values to zero.");
}

void
DataAbstract::reorderByReferenceIDs(int *reference_ids)
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
  if (m_owners.size()==2)	// Means it used to be 1 so we need to tell people
  {
	for (size_t i=0;i<m_owners.size();++i)
	{
		m_owners[i]->updateShareStatus(true);
	}
  }
/*
	for (size_t i=0;i<m_owners.size();++i)
	{
		m_owners[i]->updateShareStatus(true);
cerr << m_owners[i] << " ";
	}
cerr << endl;*/


}

void DataAbstract::removeOwner(Data* d)
{
  for (size_t i=0;i<m_owners.size();++i)
  {
	if (m_owners[i]==d)
	{
		m_owners.erase(m_owners.begin()+i,m_owners.begin()+(i+1));	// remove the element
		break;
	}
  }
  if (m_owners.size()==1)	// Means it used to be 2 so we need to tell people
  {
	m_owners[0]->updateShareStatus(false);
  }
}


}  // end of namespace
