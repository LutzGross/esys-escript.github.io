
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
#include "esysUtils/first.h"

#include "DataEmpty.h"
#include "DataException.h"


namespace {


  inline
  void
  throwStandardException(const std::string& functionName)
  {
    throw escript::DataException("Error - "+functionName+" function call invalid for DataEmpty.");
  }


  escript::DataTypes::RealVectorType dummy;	
  escript::DataTypes::CplxVectorType dummyc;	
}

namespace escript {

DataEmpty::DataEmpty() :
  parent(FunctionSpace(),DataTypes::scalarShape, true)
{

}

DataEmpty::~DataEmpty()
{
}

std::string
DataEmpty::toString() const
{
  return "(Empty Data)";
}


DataAbstract*
DataEmpty::deepCopy() const
{
  return new DataEmpty();
}

DataTypes::RealVectorType::size_type
DataEmpty::getPointOffset(int sampleNo,
                          int dataPointNo) const 
{
  throwStandardException("getPointOffset");
  return 0;
}

DataTypes::RealVectorType::size_type
DataEmpty::getPointOffset(int sampleNo,
                          int dataPointNo)
{
  throwStandardException("getPointOffset");
  return 0;
}

DataTypes::RealVectorType::size_type
DataEmpty::getLength() const
{
  return 0;
}

DataAbstract*
DataEmpty::getSlice(const DataTypes::RegionType& region) const
{
  throwStandardException("getSlice");
  return 0;
}

void
DataEmpty::setSlice(const DataAbstract* value,
                    const DataTypes::RegionType& region) 
{
  throwStandardException("setSlice");
}

int
DataEmpty::matrixInverse(DataAbstract* out) const
{
  throwStandardException("matrixInverse");
  return 0;
}


DataTypes::RealVectorType&
DataEmpty::getVectorRW()
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummy;			// dead code to stop the compiler complaining
}

const DataTypes::RealVectorType&
DataEmpty::getVectorRO() const
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummy;			// dead code to stop the compiler complaining
}


DataTypes::CplxVectorType&
DataEmpty::getVectorRWC()
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummyc;			// dead code to stop the compiler complaining
}

const DataTypes::CplxVectorType&
DataEmpty::getVectorROC() const
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummyc;			// dead code to stop the compiler complaining
}


DataTypes::RealVectorType&
DataEmpty::getTypedVectorRW(DataTypes::real_t dummypar)
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummy;			// dead code to stop the compiler complaining
}

const DataTypes::RealVectorType&
DataEmpty::getTypedVectorRO(DataTypes::real_t dummypar) const
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummy;			// dead code to stop the compiler complaining
}


DataTypes::CplxVectorType&
DataEmpty::getTypedVectorRW(DataTypes::cplx_t dummypar)
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummyc;			// dead code to stop the compiler complaining
}

const DataTypes::CplxVectorType&
DataEmpty::getTypedVectorRO(DataTypes::cplx_t dummypar) const
{
  throwStandardException("getVector");	// always throws but the compiler doesn't know that.
  return dummyc;			// dead code to stop the compiler complaining
}



void
DataEmpty::dump(const std::string fileName) const
{
    throw DataException("Error - Cannot dump() a DataEmpty object.");
}

}  // end of namespace
