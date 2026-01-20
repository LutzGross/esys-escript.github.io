
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "Data.h"
#include "WrappedArray.h"
#include "DataException.h"

#if ESYS_HAVE_NUMPY_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>
#endif

#include <iostream>

#include <boost/python/tuple.hpp>

using namespace escript;
using namespace boost::python;
using DataTypes::cplx_t;
using DataTypes::real_t;

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


// This should not be called on anything which does
// not have a __len__
bool checkForComplex(const boost::python::object& obj)
{
    try
    {
	int len=extract<int>(obj.attr("__len__")());
	for (int i=0;i<len;++i)
	{
	    const boost::python::object t=obj[i];
	    bool haslen=false;
	    try
            {
	        extract<int>(t.attr("__len__")());
	        haslen=true;
	    }
	    catch(...)
	    {
		PyErr_Clear();
	    }
	    	// If it has a length, we dig down
		// if not, we test for complex
	    if (haslen)
	    {
                if (checkForComplex(t))
		{
		    return true;
		}
	    }
	    else
	    {
        extract<DataTypes::real_t> er(t);
        if (!er.check())
		{
		    // unfortunately, if this was a numpy object, that check my fail
		    // even if it should succeed (eg numpy.int64 on python3)
		    // now we will check the dtype.kind! 
		    try
		    {
                 if ( extract<char>(t.attr("dtype").attr("kind")) == 'c' )
                {
                    return true;
                }
                else
                {
                    return false;
                }
		    }
		    catch (...)
		    {
			PyErr_Clear();
			// at this point, we have no apparent way to get a real out so 
			// we assume it must be complex
			return true;
		    }
		}
	    }
	}
	return false;
    }
    catch(...)
    {
        PyErr_Clear();
	return false;
    }
    return false;
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
:obj(obj_in),converted(false),iscomplex(false),scalar_r(nan("")),scalar_c(nan(""))
{
	dat_r=0;
	dat_c=0;
	// First we check for scalars
	try
	{
	   extract<DataTypes::cplx_t> ec(obj_in);
       extract<real_t> er(obj_in);

       if (er.check())		// check for real_t first because complex will fail this
	   {
	      scalar_r=er();
	   }
	   else
	   {
	      scalar_c=ec();
	      iscomplex=true;
	     
	   }
	   rank=0;
	   return;
	} 
	catch (...)
	{		// so we clear the failure
	   PyErr_Clear();
	}
	try
	{
	   const boost::python::object obj_in_t=obj_in[make_tuple()];
	   extract<DataTypes::cplx_t> ec(obj_in_t);
	   extract<real_t> er(obj_in_t);

	   if (er.check())
	   {	     
	      scalar_r=er();
	     
	   }
	   else
	   {
	      scalar_c=ec();
	      iscomplex=true;
	   }	   
	   rank=0;
	   return;
	} 
	catch (...)
	{		// so we clear the failure
	   PyErr_Clear();
	}


	scalar_c=0;
	scalar_r=0;
	checkFeatures(obj_in);
	getObjShape(obj,shape);
	rank=shape.size();
    iscomplex=checkForComplex(obj_in);

#if ESYS_HAVE_NUMPY_H
	// if obj is a numpy array it is much faster to copy the array through the
	// __array_struct__ interface instead of extracting single values from the
	// components via getElt(). For this to work we check below that
	// (1) this is a valid PyArrayInterface instance
	// (2) the data is stored as a contiguous C array
	// (3) the data type is suitable (correct type and byte size)
	try
	{
		object o = (extract<object>(obj.attr("__array_struct__")));
#ifdef ESPYTHON3
		if (PyCapsule_CheckExact(o.ptr()))
#else
		if (PyCObject_Check(o.ptr()))
#endif
		{
			PyObject* cobj=(PyObject*)o.ptr();
#ifdef ESPYTHON3
            const char* name = PyCapsule_GetName(cobj);
			PyArrayInterface* arr=(PyArrayInterface*)PyCapsule_GetPointer(cobj, name);
#else
			PyArrayInterface* arr=(PyArrayInterface*)PyCObject_AsVoidPtr(cobj);
#endif
#ifndef NPY_1_7_API_VERSION
  #define NPY_ARRAY_IN_ARRAY NPY_IN_ARRAY
  #define NPY_ARRAY_NOTSWAPPED NPY_NOTSWAPPED
#endif
			if (arr->two==2 && arr->flags&NPY_ARRAY_IN_ARRAY && arr->flags&NPY_ARRAY_NOTSWAPPED)
			{
				std::vector<int> strides;
				// convert #bytes to #elements
				for (int i=0; i<arr->nd; i++)
				{
					strides.push_back(arr->strides[i]/arr->itemsize);
				}

				if (arr->typekind == 'f')
				{
					if (arr->itemsize==sizeof(real_t))
					{
						convertNumpyArray<real_t>((const real_t*)arr->data, strides);
					}
					else if (arr->itemsize==sizeof(float))
			   		{
						convertNumpyArray<float>((const float*)arr->data, strides);
					}
				}
				else if (arr->typekind == 'i')
				{
					if (arr->itemsize==sizeof(int))
				   	{
						convertNumpyArray<int>((const int*)arr->data, strides);
					}
					else if (arr->itemsize==sizeof(long))
				   	{
						convertNumpyArray<long>((const long*)arr->data, strides);
					}
				}
				else if (arr->typekind == 'u')
				{
					if (arr->itemsize==sizeof(unsigned))
				   	{
						convertNumpyArray<unsigned>((const unsigned*)arr->data, strides);
					}
					else if (arr->itemsize==sizeof(unsigned long))
				   	{
						convertNumpyArray<unsigned long>((const unsigned long*)arr->data, strides);
					}
				}
				else if (arr->typekind == 'c')
				{
					if (arr->itemsize==sizeof(cplx_t))
				   	{

						convertNumpyArrayC<DataTypes::cplx_t>((const cplx_t*)arr->data, strides);
                        iscomplex=true;
					}
					// not accomodating other types of complex values
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
void WrappedArray::convertNumpyArrayC(const T* array, const std::vector<int>& strides) const
{
	// this method is only called by the constructor above which does the
	// necessary checks and initialisations
	int size=DataTypes::noValues(shape);
	dat_c=new cplx_t[size];
	switch (rank)
	{
		case 1:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				dat_c[i]=array[i*strides[0]];
			}
		break;
		case 2:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				for (int j=0;j<shape[1];j++)
				{
					dat_c[DataTypes::getRelIndex(shape,i,j)]=array[i*strides[0]+j*strides[1]];
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
						dat_c[DataTypes::getRelIndex(shape,i,j,k)]=array[i*strides[0]+j*strides[1]+k*strides[2]];
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
							dat_c[DataTypes::getRelIndex(shape,i,j,k,m)]=array[i*strides[0]+j*strides[1]+k*strides[2]+m*strides[3]];
						}
					}
				}
			}
		break;
	}
}


template<typename T>
void WrappedArray::convertNumpyArray(const T* array, const std::vector<int>& strides) const
{
	// this method is only called by the constructor above which does the
	// necessary checks and initialisations
	int size=DataTypes::noValues(shape);
	dat_r=new real_t[size];
	switch (rank)
	{
		case 1:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				dat_r[i]=array[i*strides[0]];
			}
		break;
		case 2:
#pragma omp parallel for
			for (int i=0;i<shape[0];i++)
			{
				for (int j=0;j<shape[1];j++)
				{
					dat_r[DataTypes::getRelIndex(shape,i,j)]=array[i*strides[0]+j*strides[1]];
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
						dat_r[DataTypes::getRelIndex(shape,i,j,k)]=array[i*strides[0]+j*strides[1]+k*strides[2]];
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
							dat_r[DataTypes::getRelIndex(shape,i,j,k,m)]=array[i*strides[0]+j*strides[1]+k*strides[2]+m*strides[3]];
						}
					}
				}
			}
		break;
	}
}


void WrappedArray::convertArrayR() const
{
	if ((converted) || (rank<=0) || (rank>4))	// checking illegal rank here to avoid memory issues later
	{					// yes the failure is silent here but not doing the copy 
	    return;				// will just cause an error to be raised later
	}
	int size=DataTypes::noValues(shape);
	real_t* tdat=new real_t[size];
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
	dat_r=tdat;    
	converted=true;
}  


void WrappedArray::convertArrayC() const
{
	if ((converted) || (rank<=0) || (rank>4))	// checking illegal rank here to avoid memory issues later
	{					// yes the failure is silent here but not doing the copy 
	    return;				// will just cause an error to be raised later
	}
	int size=DataTypes::noValues(shape);
	cplx_t* tdat=new cplx_t[size];
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
	dat_c=tdat;    
	converted=true;
}  


void WrappedArray::convertArray() const
{
    if (iscomplex)
    {
	convertArrayC();
    }
    else
    {
	convertArrayR();
    }
}

WrappedArray::~WrappedArray()
{
    if (dat_r!=0)
    {
	delete[] dat_r;
    }
    if (dat_c!=0)
    {
	delete[] dat_c;
    }
}


