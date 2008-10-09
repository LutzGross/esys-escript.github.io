
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


#include "DataLazy.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif
#include "FunctionSpace.h"
#include "DataTypes.h"

using namespace std;

namespace escript
{

const std::string&
opToString(ES_optype op);

namespace
{

// return the FunctionSpace of the result of "left op right"
FunctionSpace
resultFS(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
	// perhaps this should call interpolate and throw or something?
	// maybe we need an interpolate node -
	// that way, if interpolate is required in any other op we can just throw a 
	// programming error exception.
	return left->getFunctionSpace();
}

// return the shape of the result of "left op right"
DataTypes::ShapeType
resultShape(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
	return DataTypes::scalarShape;
}

size_t
resultLength(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
{
   switch(op)
   {
   case IDENTITY: return left->getLength();
   default: 
	throw DataException("Programmer Error - attempt to getLength() for operator "+opToString(op)+".");

   }
}

string ES_opstrings[]={"UNKNOWN","IDENTITY"};
int ES_opcount=2;

}	// end anonymous namespace


const std::string&
opToString(ES_optype op)
{
  if (op<0 || op>=ES_opcount) 
  {
    op=UNKNOWNOP;
  }
  return ES_opstrings[op];
}


DataLazy::DataLazy(DataAbstract_ptr p)
	: parent(p->getFunctionSpace(),p->getShape()),
	m_left(p),
	m_op(IDENTITY)
{
   length=resultLength(m_left,m_right,m_op);
}

DataLazy::DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op)
	: parent(resultFS(left,right,op), resultShape(left,right,op)),
	m_left(left),
	m_right(right),
	m_op(op)
{
   length=resultLength(m_left,m_right,m_op);
}


DataLazy::~DataLazy()
{
}

// If resolving records a pointer to the resolved Data we may need to rethink the const on this method
DataReady_ptr
DataLazy::resolve()
{
  if (m_op==IDENTITY)
  {
    return m_left->resolve();
  }
  else
  {
    throw DataException("Programmer error - do not know how to resolve operator "+opToString(m_op)+".");
  }
}

std::string
DataLazy::toString() const
{
  return "Lazy evaluation object. No details available.";
}

// Note that in this case, deepCopy does not make copies of the leaves.
// Hopefully copy on write (or whatever we end up using) will take care of this.
DataAbstract* 
DataLazy::deepCopy()
{
  if (m_op==IDENTITY)
  {
	return new DataLazy(m_left);	// we don't need to copy the child here
  }
  return new DataLazy(m_left->deepCopy()->getPtr(),m_right->deepCopy()->getPtr(),m_op); 
}


DataTypes::ValueType::size_type
DataLazy::getLength() const
{
  return length;
}


DataAbstract*
DataLazy::getSlice(const DataTypes::RegionType& region) const
{
  // this seems like a really good one to include I just haven't added it yet
  throw DataException("getSlice - not implemented for Lazy objects - yet.");
}

DataTypes::ValueType::size_type 
DataLazy::getPointOffset(int sampleNo,
                 int dataPointNo) const
{
  throw DataException("getPointOffset - not implemented for Lazy objects - yet.");
}

}	// end namespace
