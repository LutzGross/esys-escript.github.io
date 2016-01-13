
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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

#include <iostream>

using namespace escript;
using namespace boost::python;

namespace
{

void checkFeatures(const boost::python::object& obj)
{
	using namespace std;
	int len=0;
	boost::python::object o2;
	try
	{
	   len=extract<int>(obj.attr("__len__")());
	}
	catch (...)
	{
	   PyErr_Clear();
	   string s=extract<std::string>(obj.attr("__str__")());
	   cerr << "Failure for object: " << s << endl;
	   s=extract<std::string>(obj.attr("__repr__")());
	   cerr << "Failurr for object: " << s << endl;
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
	// First we check for scalars
	try
	{
	   double v=extract<double>(obj_in);
	   m_scalar=v;
	   rank=0;
	   return;
	} 
	catch (...)
	{		// so we 
	   PyErr_Clear();
	}
	try
	{
	   double v=extract<double>(obj_in[make_tuple()]);
	   m_scalar=v;
	   rank=0;
	   return;
	} 
	catch (...)
	{		// so we 
	   PyErr_Clear();
	}


	m_scalar=0;
	checkFeatures(obj_in);
	getObjShape(obj,shape);
	rank=shape.size();
}



