// $Id$
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

#include "DataFactory.h"

#include <boost/python/extract.hpp>

using namespace boost::python;

namespace escript {

Data
Scalar(double value,
       const FunctionSpace& what,
       bool expanded)
{
    //
    // an empty shape is a scalar
    DataArrayView::ShapeType shape;
    return Data(value,shape,what,expanded);
}

Data
Vector(double value,
       const FunctionSpace& what,
       bool expanded)
{
    DataArrayView::ShapeType shape(1,what.getDomain().getDim());
    return Data(value,shape,what,expanded);
}

Data
Tensor(double value,
       const FunctionSpace& what,
       bool expanded)
{
    DataArrayView::ShapeType shape(2,what.getDomain().getDim());
    return Data(value,shape,what,expanded);
}

Data
Tensor3(double value,
        const FunctionSpace& what,
        bool expanded)
{
    DataArrayView::ShapeType shape(3,what.getDomain().getDim());
    return Data(value,shape,what,expanded);
}

Data
Tensor4(double value,
        const FunctionSpace& what,
        bool expanded)
{
    DataArrayView::ShapeType shape(4,what.getDomain().getDim());
    return Data(value,shape,what,expanded);
}

Data
convertToData(const boost::python::object& value,
              const FunctionSpace& what) 
{
     // first we try to extract a Data object from value 
     extract<Data> value_data(value);
     if (value_data.check()) {
         Data extracted_data=value_data();
         if (extracted_data.isEmpty()) {
            return extracted_data;
         } else {
            return Data(extracted_data,what);
         }
     } else {
        return Data(value,what);
     }
}

}  // end of namespace
