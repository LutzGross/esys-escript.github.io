
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined escript_DataFactory_20040721_H
#define escript_DataFactory_20040721_H

#include "system_dep.h"

#include "AbstractDomain.h"
#include "FunctionSpace.h"
#include "Data.h"

#include <boost/python/object.hpp>

namespace escript {

/**
   \brief
    A collection of factory functions for creating Data objects which contain
    data points of various shapes.
*/

/**
   \brief
   Return a Data object containing scalar data-points.
   ie: rank 0 data-points.
   \param value - Input - Single value applied to all Data.
   \param what - Input - A description of what this data represents.
   \param expanded - Input - if true fill the entire container with 
                     the value. Otherwise a more efficient storage 
                     mechanism will be used.
*/
ESCRIPT_DLL_API Data
Scalar(double value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);


/**
   \brief
   Return a Data object containing scalar data-points.
   ie: rank 0 data-points.
   \param value - Input - Single value applied to all Data.
   \param what - Input - A description of what this data represents.
   \param expanded - Input - if true fill the entire container with 
                     the value. Otherwise a more efficient storage 
                     mechanism will be used.
*/
ESCRIPT_DLL_API Data
Scalar(DataTypes::cplx_t value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

/**
 * \brief Python interface for the Scalar factory 
*/
ESCRIPT_DLL_API
Data
ScalarFromObj(boost::python::object o,
	const FunctionSpace& what=FunctionSpace(),
	bool expanded=false);


/**
   \brief
   Return a Data object containing vector data-points.
   ie: rank 1 data-points.
*/
ESCRIPT_DLL_API Data
Vector(double value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

ESCRIPT_DLL_API
Data
VectorFromObj(boost::python::object o,
	const FunctionSpace& what=FunctionSpace(),
	bool expanded=false);

/**
   \brief
   Return a Data object containing tensor datapoints.
   ie: rank 2 data-points.
*/
ESCRIPT_DLL_API Data
Tensor(double value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

ESCRIPT_DLL_API
Data
TensorFromObj(boost::python::object o,
	const FunctionSpace& what=FunctionSpace(),
	bool expanded=false);
/**
   \brief
   Return a Data object containing tensor3 datapoints.
   ie: rank 3 data-points.
*/
ESCRIPT_DLL_API Data
Tensor3(double value,
        const FunctionSpace& what=FunctionSpace(),
        bool expanded=false);

ESCRIPT_DLL_API
Data
Tensor3FromObj(boost::python::object o,
	const FunctionSpace& what=FunctionSpace(),
	bool expanded=false);

/**
   \brief
   Return a Data object containing tensor4 datapoints.
   ie: rank 4 data-points.
*/
ESCRIPT_DLL_API Data
Tensor4(double value,
        const FunctionSpace& what=FunctionSpace(),
        bool expanded=false);

ESCRIPT_DLL_API
Data
Tensor4FromObj(boost::python::object o,
	const FunctionSpace& what=FunctionSpace(),
	bool expanded=false);

/**
   \brief
   reads Data on domain from file in netCDF format
*/
ESCRIPT_DLL_API Data 
load(const std::string fileName,
     const AbstractDomain& domain);
/**
   \brief
   returns true if the load funtion is configured.
*/
ESCRIPT_DLL_API bool
loadConfigured();

/**
   \brief
   Tries to convert value into a Data object on FunctionSpace what.
   If value is already a Data object, the object is returned if it is defined on what otherwise
   interpolated data of values are returned. If value is not a data object it is tried to generate
   the corresponding data object. escript::DataEmpty() is returned if value is identified as empty.
*/
ESCRIPT_DLL_API Data
convertToData(const boost::python::object& value,
              const FunctionSpace& what=FunctionSpace());


} // end of namespace

#endif
