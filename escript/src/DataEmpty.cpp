
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

#include "DataEmpty.h"
#include "DataException.h"

namespace escript {

DataEmpty::DataEmpty() :
  DataAbstract(FunctionSpace(),DataTypes::scalarShape)
{
//  resetPointDataView();
}

DataEmpty::~DataEmpty()
{
}

std::string
DataEmpty::toString() const
{
  return "(Empty Data)";
}

DataTypes::ValueType::size_type
DataEmpty::getPointOffset(int sampleNo,
                          int dataPointNo) const 
{
  throwStandardException("getPointOffset");
  return 0;
}

// DataArrayView
// DataEmpty::getDataPoint(int sampleNo,
//                         int dataPointNo)
// {
//   throwStandardException("getDataPoint");
//   return getPointDataView();
// }

DataTypes::ValueType::size_type
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

void
DataEmpty::throwStandardException(const std::string& functionName) const
{
  throw DataException("Error - "+functionName+" function call invalid for DataEmpty.");
}

DataTypes::ValueType&
DataEmpty::getVector()
{
  throwStandardException("getVector");
}

const DataTypes::ValueType&
DataEmpty::getVector() const
{
  throwStandardException("getVector");
}



}  // end of namespace
