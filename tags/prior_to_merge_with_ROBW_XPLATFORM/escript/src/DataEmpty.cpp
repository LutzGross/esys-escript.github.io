//$Id$

/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#include "DataEmpty.h"
#include "DataException.h"

namespace escript {

DataEmpty::DataEmpty() :
  DataAbstract(FunctionSpace())
{
  resetPointDataView();
}

DataEmpty::~DataEmpty()
{
}

std::string
DataEmpty::toString() const
{
  return "(Empty Data)";
}

DataArrayView::ValueType::size_type
DataEmpty::getPointOffset(int sampleNo,
                          int dataPointNo) const 
{
  throwStandardException("getPointOffset");
  return 0;
}

DataArrayView
DataEmpty::getDataPoint(int sampleNo,
                        int dataPointNo)
{
  throwStandardException("getDataPoint");
  return getPointDataView();
}

DataArrayView::ValueType::size_type
DataEmpty::getLength() const
{
  return 0;
}

DataAbstract*
DataEmpty::getSlice(const DataArrayView::RegionType& region) const
{
  throwStandardException("getSlice");
  return 0;
}

void
DataEmpty::setSlice(const DataAbstract* value,
                    const DataArrayView::RegionType& region) 
{
  throwStandardException("setSlice");
}

void
DataEmpty::reshapeDataPoint(const DataArrayView::ShapeType& shape)
{
  throwStandardException("reshapeDataPoint");
}

void
DataEmpty::throwStandardException(const std::string& functionName) const
{
  throw DataException("Error - "+functionName+" function call invalid for DataEmpty.");
}

}  // end of namespace
