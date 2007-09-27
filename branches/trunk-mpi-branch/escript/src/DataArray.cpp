
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

#include "DataArray.h"

#include <boost/python/extract.hpp>

using namespace boost::python;
using namespace std;

namespace escript {

DataArray::DataArray(double value)
{
    m_data.resize(1,value,1);
    // create a view with an empty shape, a scalar.
    m_dataView.reset(new DataArrayView(m_data,DataArrayView::ShapeType()));
}

DataArray::DataArray(const DataArrayView::ShapeType& shape,
                     double value)
{
    int len = DataArrayView::noValues(shape);
    m_data.resize(len,value,len);
    m_dataView.reset(new DataArrayView(m_data,shape));
}

DataArray::DataArray(const DataArray& value):
    m_data(value.getData())
{
    m_dataView.reset(new DataArrayView(m_data,value.getView().getShape()));
}

DataArray::DataArray(const DataArrayView& value):
    m_data(value.getData())
{
    m_dataView.reset(new DataArrayView(m_data,value.getShape()));
}

DataArray::DataArray(const object& value)
{
    // this will throw if the value cannot be represented
    numeric::array asNumArray(value);
    initialise(asNumArray);
}

DataArray::DataArray(const boost::python::numeric::array& value) 
{
    initialise(value);
}
 
void
DataArray::initialise(const boost::python::numeric::array& value)
{
    // extract the shape of the numarray
    DataArrayView::ShapeType tempShape;  
    for (int i=0; i<value.getrank(); i++) {
      tempShape.push_back(extract<int>(value.getshape()[i]));
    }
    // allocate the space for the data vector
    int len = DataArrayView::noValues(tempShape);
    m_data.resize(len,0.,len);
    // create a view with the same shape
    m_dataView.reset(new DataArrayView(m_data,tempShape));
    // fill the data vector with the values from the numarray
    m_dataView->copy(value);
}

const DataArrayView&
DataArray::getView() const 
{
    return *m_dataView;
}

DataArrayView&
DataArray::getView()
{
    return *m_dataView;
}

const DataArrayView::ValueType&
DataArray::getData() const
{
    return m_data;
}

DataArrayView::ValueType&
DataArray::getData()
{
    return m_data;
}

}  // end of namespace
