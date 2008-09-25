
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

using namespace std;

namespace escript {

DataAbstract::DataAbstract(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty):
    m_noDataPointsPerSample(what.getNumDPPSample()),
    m_noSamples(what.getNumSamples()),
    m_functionSpace(what),
    m_shape(shape),
    m_rank(DataTypes::getRank(shape)),
    m_novalues(DataTypes::noValues(shape))
{
    m_isempty=isDataEmpty;
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
    if (!((right.getRank()==0) || 
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
DataAbstract::copyAll(const boost::python::numeric::array& value)
{
    throw DataException("Error - DataAbstract::copying data from numarray objects is not supported.");
}
void
DataAbstract::copyToDataPoint(const int sampleNo, const int dataPointNo, const double value)
{
    throw DataException("Error - DataAbstract::copying data from double value to a single data point is not supported.");
}
void
DataAbstract::copyToDataPoint(const int sampleNo, const int dataPointNo, const boost::python::numeric::array& value)
{
    throw DataException("Error - DataAbstract::copying data from numarray objects to a single data point is not supported.");
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


DataTypes::ValueType&
DataAbstract::getVector()
{
   throw DataException("Error - DataAbstract:: does not have a DataVector.");
}

const DataTypes::ValueType&
DataAbstract::getVector() const
{
   throw DataException("Error - DataAbstract:: does not have a DataVector.");
}




}  // end of namespace
