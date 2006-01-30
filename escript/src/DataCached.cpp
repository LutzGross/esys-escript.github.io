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

#include "DataCached.h"
#include "DataException.h"

namespace escript {

DataCached::DataCached() :
  DataAbstract(FunctionSpace())
{
  resetPointDataView();
}

DataCached::~DataCached()
{
}

std::string
DataCached::toString() const
{
  return "(Cached Data)";
}

DataArrayView::ValueType::size_type
DataCached::getPointOffset(int sampleNo,
                          int dataPointNo) const 
{
  throwStandardException("getPointOffset");
  return 0;
}

DataArrayView
DataCached::getDataPoint(int sampleNo,
                        int dataPointNo)
{
  throwStandardException("getDataPoint");
  return getPointDataView();
}

DataArrayView::ValueType::size_type
DataCached::getLength() const
{
  return 0;
}

DataAbstract*
DataCached::getSlice(const DataArrayView::RegionType& region) const
{
  throwStandardException("getSlice");
  return 0;
}

void
DataCached::setSlice(const DataAbstract* value,
                    const DataArrayView::RegionType& region) 
{
  throwStandardException("setSlice");
}

void
DataCached::reshapeDataPoint(const DataArrayView::ShapeType& shape)
{
  throwStandardException("reshapeDataPoint");
}

void
DataCached::throwStandardException(const std::string& functionName) const
{
  throw DataException("Error - "+functionName+" function call invalid for DataCached.");
}

}  // end of namespace
