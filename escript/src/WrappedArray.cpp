
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

#include "WrappedArray.h"
#include "DataException.h"

using namespace escript;
using namespace boost::python;

namespace
{

void checkFeatures(const boost::python::object& obj)
{
	int len=0;
	boost::python::object o2;
	try
	{
	   len=extract<int>(obj.attr("__len__")());
	}
	catch (...)
	{
	   PyErr_Clear();
	   throw DataException("Object passed to WrappedArray must support __len__");
	}
	try
	{
	   o2=obj.attr("__getitem__");
	}
	catch (...)
	{
	   PyErr_Clear();
	   throw DataException("Object passed to WrappedArray must support __getitem__");
	}
}

void getObjShape(const boost::python::object& obj, DataTypes::ShapeType& s)
{
	int len=0;
	try
	{
	   len=extract<int>(obj.attr("__len__")());
	}
	catch(...)
	{
	   PyErr_Clear();		// tell python the error isn't there anymore
	   return;
	}
	if (len<1)
	{
	   throw DataException("Array filter - no empty components in arrays please.");
	}
	s.push_back(len);	

	if (s.size()>ESCRIPT_MAX_DATA_RANK)
	{
	   throw DataException("Array filter - Maximum rank exceeded in array");
	}
	getObjShape(obj[0],s);
}

}

WrappedArray::WrappedArray(const boost::python::object& obj_in)
:obj(obj_in)
{
	checkFeatures(obj_in);
	getObjShape(obj,shape);
	rank=shape.size();
}



