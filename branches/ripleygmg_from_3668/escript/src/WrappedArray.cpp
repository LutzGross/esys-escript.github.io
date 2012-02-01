
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <boost/python/tuple.hpp>
#include "WrappedArray.h"
#include "DataException.h"
#if HAVE_NUMPY_H
#include <numpy/ndarrayobject.h>
#endif

#include <iostream>

using namespace escript;
using namespace boost::python;

namespace
{

void checkFeatures(const boost::python::object& obj)
{
	using namespace std;
	boost::python::object o2;
	try
	{
	   /*int len=*/ extract<int>(obj.attr("__len__")());
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
	dat=0;
	// First we check for scalars
	try
	{
	   double v=extract<double>(obj_in);
	   m_scalar=v;
	   rank=0;
	   return;
	} 
	catch (...)
	{		// so we clear the failure
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
	{		// so we clear the failure
	   PyErr_Clear();
	}


	m_scalar=0;
	checkFeatures(obj_in);
	getObjShape(obj,shape);
	rank=shape.size();

#if HAVE_NUMPY_H
	// if obj is a numpy array it is much faster to copy the array through the
	// __array_struct__ interface instead of extracting single values from the
	// components via getElt(). For this to work we check below that
	// (1) this is a valid PyArrayInterface instance
	// (2) the data is stored as a contiguous C array
	// (3) the data type is suitable (correct type and byte size)
	try
	{
		object o = (extract<object>(obj.attr("__array_struct__")));
		if (PyCObject_Check(o.ptr()))
		{
			PyObject* cobj=(PyObject*)o.ptr();
			PyArrayInterface* arr=(PyArrayInterface*)PyCObject_AsVoidPtr(cobj);
			if (arr->two==2 && arr->flags&NPY_IN_ARRAY)
			{
				if (arr->typekind == 'f' && arr->itemsize==sizeof(double))
			   	{
					convertNumpyArray<double>((const double*)arr->data);
				}
				else if (arr->typekind == 'f' && arr->itemsize==sizeof(float))
			   	{
					convertNumpyArray<float>((const float*)arr->data);
				}
				else if (arr->typekind == 'i' && arr->itemsize==sizeof(int))
			   	{
					convertNumpyArray<int>((const int*)arr->data);
				}
			}
		}
	} catch (...)
	{
		PyErr_Clear();
	}
#endif
}

template<typename T>
void WrappedArray::convertNumpyArray(const T* array) const
{
	// this method is only called by the constructor above which does the
	// necessary checks and initialisations
	int size=DataTypes::noValues(shape);
	dat=new double[size];
	switch (rank)
	{
		case 1:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				dat[i]=array[i];
			}
		break;
		case 2:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				for (int j=0;j<shape[1];j++)
				{
					dat[DataTypes::getRelIndex(shape,i,j)]=array[j+i*shape[1]];
				}
			}
		break;
		case 3:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				for (int j=0;j<shape[1];j++)
				{
					for (int k=0;k<shape[2];k++)
					{
						dat[DataTypes::getRelIndex(shape,i,j,k)]=array[k+j*shape[2]+i*shape[1]*shape[2]];
					}
				}
			}
		break;
		case 4:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				for (int j=0;j<shape[1];j++)
				{
					for (int k=0;k<shape[2];k++)
					{
						for (int m=0;m<shape[3];m++)
						{
							dat[DataTypes::getRelIndex(shape,i,j,k,m)]=array[m+k*shape[3]+j*shape[2]*shape[3]+i*shape[1]*shape[2]*shape[3]];
						}
					}
				}
			}
		break;
	}
}

void WrappedArray::convertArray() const
{
	if ((dat!=0) || (rank<=0) || (rank>4))	// checking illegal rank here to avoid memory issues later
	{					// yes the failure is silent here but not doing the copy 
	    return;				// will just cause an error to be raised later
	}
	int size=DataTypes::noValues(shape);
	double* tdat=new double[size];
	switch (rank)
	{
	case 1: for (int i=0;i<shape[0];i++)
		{
			tdat[i]=getElt(i);
		}
		break;
	case 2: for (int i=0;i<shape[0];i++)
		{
		    for (int j=0;j<shape[1];j++)
		    {
			tdat[DataTypes::getRelIndex(shape,i,j)]=getElt(i,j);
		    }
		}
		break;
	case 3: for (int i=0;i<shape[0];i++)
		{
		    for (int j=0;j<shape[1];j++)
		    {
			for (int k=0;k<shape[2];k++)
			{
			    tdat[DataTypes::getRelIndex(shape,i,j,k)]=getElt(i,j,k);
			}
		    }
		}
		break;
	case 4: for (int i=0;i<shape[0];i++)
		{
		    for (int j=0;j<shape[1];j++)
		    {
			for (int k=0;k<shape[2];k++)
			{
			    for (int m=0;m<shape[3];m++)
			    {
			    	tdat[DataTypes::getRelIndex(shape,i,j,k,m)]=getElt(i,j,k,m);
			    }
			}
		    }
		}
		break;
	default:
		;  // do nothing
		// can't happen. We've already checked the bounds above
	}
	dat=tdat;
}

WrappedArray::~WrappedArray()
{
	if (dat!=0)
	{
	    delete dat;
	}
}


