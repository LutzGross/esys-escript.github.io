//$Id$
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "escript/Data/DataEmpty.h"
#include "escript/Data/DataException.h"

namespace escript {

DataEmpty::DataEmpty() :
  DataAbstract(FunctionSpace())
{
}

DataEmpty::~DataEmpty()
{
}

std::string
DataEmpty::toString() const
{
  return "Empty Data.";
}

void
DataEmpty::throwStandardException(const std::string& functionName) const
{
  throw DataException("Error - "+functionName+" function call invalid for DataEmpty.");
}

DataArrayView::ValueType::size_type
DataEmpty::getPointOffset(int sampleNo,
                          int dataPointNo) const 
{
  throwStandardException("getPointOffset");
  return 0;
}

double*
DataEmpty::getSampleDataByTag(int tag)
{
  throwStandardException("getSampleDataByTag");
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

}  // end of namespace
